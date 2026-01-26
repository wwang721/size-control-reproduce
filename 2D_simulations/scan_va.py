import numpy as np
from tqdm import tqdm
from mpi4py import MPI

# for computing connected components
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def compute_forces(positions, connect, L, rc, lb):
    '''Function to compute spring forces (vectorized)'''
    N = positions.shape[0]  # Number of particles

    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=lb)  # cutoff length
    pairs_array = np.array(list(pairs), dtype=np.intp).reshape(-1, 2)     # shape (M, 2)

    adj_p = np.zeros((N, N), dtype=bool)
    adj_c = np.zeros((N, N), dtype=bool)

    # fill adjacency
    adj_p[pairs_array[:,0], pairs_array[:,1]] = True
    adj_c[connect[:,0], connect[:,1]] = True

    # now a single fancy‐index lookup gives you the masks:
    mask_c = adj_p[connect[:,0],    connect[:,1]]
    mask_p = adj_c[pairs_array[:,0], pairs_array[:,1]]

    # Pairs in connect but not in current neighbor list
    mask_missing = ~mask_c
    num_breaking = np.count_nonzero(mask_missing)

    if pairs_array.shape[0] == 0:
        # No pairs found, return zero forces
        return np.zeros_like(positions), pairs_array, num_breaking


    i_idx = pairs_array[:, 0]
    j_idx = pairs_array[:, 1]

    # Compute displacements between neighbors
    # (M, 2) where M is the number of pairs
    rij = positions[j_idx] - positions[i_idx]
    dist = np.linalg.norm(rij, axis=1)

    # Mask for interacting pairs
    mask = dist < rc  # rc is the contact distance
    
    # Include all historical connections in mask
    mask = np.logical_or(mask, mask_p)

    i_masked = i_idx[mask]
    j_masked = j_idx[mask]
    rij = rij[mask]
    dist = dist[mask]

    # Compute spring pair forces
    force_mags = dist - L
    force_dirs = rij / dist[:, np.newaxis]  # (M, 2)
    forces = force_mags[:, np.newaxis] * force_dirs  # (M, 2)

    # Accumulate forces
    net_forces = np.zeros_like(positions)
    np.add.at(net_forces, i_masked, forces)
    np.add.at(net_forces, j_masked, -forces)

    return net_forces, pairs_array[mask], num_breaking


def rebuild_connection_matrix(N, connect):
    i, j = connect[:, 0], connect[:, 1]
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    data = np.ones(len(rows), dtype=bool)
    connect_matrix = csr_matrix((data, (rows, cols)), shape=(N, N))

    return connect_matrix


def cell_division(dividing_cells, positions, va, va_x, va_y, connect, daughter_distance):
    N = positions.shape[0]  # number of particles

    connect_matrix = rebuild_connection_matrix(N, connect)

    # Coordination number: degree = row sum
    # coordinations = np.diff(connect_matrix.indptr)

    for idx_dividing in dividing_cells:
        # Get the indices of the connected cells
        connected_cells = connect_matrix.indices[connect_matrix.indptr[idx_dividing]: connect_matrix.indptr[idx_dividing+1]]

        angle_dividing = np.random.uniform(-np.pi/2., np.pi/2.)
        direction = np.array([np.cos(angle_dividing), np.sin(angle_dividing)])

        pos_divide = positions[idx_dividing]

        offset = (daughter_distance/2.) * direction
        pos1 = pos_divide - offset
        pos2 = pos_divide + offset

        # Positions of connected cells
        connected_pos = positions[connected_cells]

        # Compute distances to each daughter
        dist1 = np.linalg.norm(connected_pos - pos1, axis=1)
        dist2 = np.linalg.norm(connected_pos - pos2, axis=1)

        # Decide assignments
        assign_to_1 = dist1 < dist2
        daughter1_neighbors = connected_cells[assign_to_1]
        daughter2_neighbors = connected_cells[~assign_to_1]

        # Update positions
        positions[idx_dividing] = pos1  # daughter1
        positions = np.concatenate([positions, [pos2]], axis=0)  # daughter2

        # Update active velocities
        va_x[idx_dividing] = va * np.cos(angle_dividing - np.pi)
        va_y[idx_dividing] = va * np.sin(angle_dividing - np.pi)

        va_x = np.concatenate([va_x, [va * np.cos(angle_dividing)]])
        va_y = np.concatenate([va_y, [va * np.sin(angle_dividing)]])
        

        # Update connectivity
        # Convert to LIL format for efficient row/col updates
        connect_matrix = connect_matrix.tolil()

        # Resize to (N+1, N+1)
        connect_matrix.resize((N + 1, N + 1))

        # Remove old connections
        connect_matrix[idx_dividing, connected_cells] = 0
        connect_matrix[connected_cells, idx_dividing] = 0

        # Add new connections
        connect_matrix[idx_dividing, daughter1_neighbors] = 1
        connect_matrix[daughter1_neighbors, idx_dividing] = 1
        connect_matrix[N, daughter2_neighbors] = 1
        connect_matrix[daughter2_neighbors, N] = 1

        # Convert back to CSR if needed
        connect_matrix = connect_matrix.tocsr()

        N += 1

    row, col = connect_matrix.nonzero()
    mask = row < col
    connect = np.stack([row[mask], col[mask]], axis=1)

    return positions, va_x, va_y, connect


def select_daughter_cluster(N, connect):
    """Select a daughter cluster randomly from the connected components."""
    n_components, labels = connected_components(csgraph=rebuild_connection_matrix(N, connect), directed=False)

    # Randomly select a daughter cluster
    if n_components > 1:
        selected_idx = np.random.randint(n_components)
        selected_cells = np.where(labels == selected_idx)[0]
        N = len(selected_cells)

        selected_mask = np.isin(connect[:, 0], selected_cells) & np.isin(connect[:, 1], selected_cells)
        selected_connect = connect[selected_mask]

        # Use searchsorted to map global → local index
        selected_connect_local = np.searchsorted(selected_cells, selected_connect)

        return selected_cells, N, selected_connect_local
    else:
        return None, N, connect


# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

vas = np.arange(0.55, 1.51, 0.01)[::-1]   # Self-propulsion velocity - 30 microns/h
# vas = np.linspace(0.5, 1, 51)[::-1]    
# vas = np.linspace(1.01, 1.5, 50)


counter = 0
for va in vas:
    if rank == 0:
        print(
            f"\nRunning simulation for va={va}: {counter}/{len(vas)}")

    counter += 1

    M = 1*100        # (# of local tasks * size) Total number of realizations

    # Parameters
    # time unit: 1 hour
    # length unit: 40 microns
    dt = 0.01        # Time step
    T = 1000       # Total simulation time
    steps = int(T / dt)     # Number of steps
    mu = 1.0        # Mobility
    k = 4.0     # Spring constant - 1/\mu k = 15 min
    L = 1.0     # relaxation length - 40 microns

    # Dr = 1.33     # Noise strength - 1/45min
    tau = 0.5    # AOUP relaxation time - 0.5 hr
    # va = 0.5   # Self-propulsion velocity - 30 microns/h - v_rms for AOUPs
    D = (va**2) * tau / 2. # AOUP diffusion coefficient - in 2D: v_rms = sqrt(2*D/tau)

    kd = np.log(2)/16.     # Cell division rate - 16-hour cell cycle
    rc = 1.001*L  # Contact distance
    lb = 1.1625*L  # Breaking length - stretch 6.5 microns



    # Split realizations among processes
    M_local = M // size  # Number of realizations per process
    if rank < M % size:
        M_local += 1

    if rank == 0:
        pbar = tqdm(total=M_local, desc="Processing")

    update_interval = steps // 20000  # update 10000 times total

    N_local = []
    breaking_local = []
    connections_local = []
    rupture_times_local = []

    # Process realizations
    for m in range(M_local):

        N = 1        # Number of particles
        N_trajectory = [N]  # Store the number of particles at each step
        t_trajectory = [0.]  # Store the time at each step
        rupture_times = [0.]  # Store rupture times

        # Initialize positions and velocities
        xx = (np.arange(np.sqrt(N)) - (np.sqrt(N)-1)/2.)*L
        x, y = np.meshgrid(xx, xx)
        x = x.flatten()
        y = y.flatten()
        positions = np.array([x, y]).T
        
        connect = np.array([], dtype=np.intp).reshape(-1, 2)

        # Initialize the active velocities
        va_x = np.sqrt(D/tau) * np.random.randn(N)
        va_y = np.sqrt(D/tau) * np.random.randn(N)

        tot_num_breaking = 0
        tot_num_connections = 0

        # Simulation loop
        for t in range(steps):
            forces, connect, num_breaking = compute_forces(positions, connect, L, rc, lb)

            num_connections = connect.shape[0] # number of connections
            tot_num_connections += num_connections
            tot_num_breaking += num_breaking

            # randomly select a daughter cluster
            selected, N, connect = select_daughter_cluster(N, connect)
            if selected is not None:
                positions = positions[selected]
                x = positions[:, 0]
                y = positions[:, 1]
                va_x = va_x[selected]
                va_y = va_y[selected]
                forces = forces[selected]
                
                rupture_times.append((t+1)*dt)

            # Update positions and angles using Euler's method
            vx = mu * k * forces[:, 0] + va_x
            vy = mu * k * forces[:, 1] + va_y

            # update position first to ensure Euler's method
            x += vx * dt
            y += vy * dt
            positions = np.array([x, y]).T

            d_va_x = - va_x * dt / tau
            d_va_y = - va_y * dt / tau
            # Gaussian white noise
            noise_x = np.sqrt(2 * D * dt) * np.random.randn(N) / tau
            noise_y = np.sqrt(2 * D * dt) * np.random.randn(N) / tau
            va_x += d_va_x + noise_x
            va_y += d_va_y + noise_y

            # cell division
            dividing_cells = np.nonzero(np.random.rand(N) < kd * dt)[0]
            if len(dividing_cells) > 0:
                positions, va_x, va_y, connect = cell_division(dividing_cells, positions, va, va_x, va_y, connect, 2*L/3.)
                x = positions[:, 0]
                y = positions[:, 1]
                N = len(x)

            if t % update_interval == 0:
                N_trajectory.append(N)
                t_trajectory.append((t+1)*dt)

        if rank == 0:
            pbar.update(1)

        N_local.append(N_trajectory)
        breaking_local.append(tot_num_breaking)
        connections_local.append(tot_num_connections*dt)
        rupture_times_local.append(rupture_times)


    # Gather results from all processes
    N_all = comm.gather(N_local, root=0)
    breaking_all = comm.gather(breaking_local, root=0)
    connections_all = comm.gather(connections_local, root=0)
    rupture_times_all = comm.gather(rupture_times_local, root=0)

    # Combine results on the root process
    if rank == 0:
        N_trajectory_all = np.concatenate(N_all)
        tot_num_breaking_all = np.concatenate(breaking_all)
        tot_num_connection_hours_all = np.concatenate(connections_all)
        rupture_times_all = np.concatenate([r_times for r_process in rupture_times_all for r_times in r_process])
        np.savez('data/va%g.npz' % va, t=t_trajectory, N=N_trajectory_all, breaking=tot_num_breaking_all, connections=tot_num_connection_hours_all, rupture=rupture_times_all)
        
    if rank == 0:
        pbar.close()
        print('\nSaved data to data/va%g.npz' % va)
        
if rank == 0:
    print("Done.")
