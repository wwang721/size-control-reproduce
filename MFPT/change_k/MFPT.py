import numpy as np
from tqdm import tqdm
from mpi4py import MPI

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Parameters
M = 1000             # Number of realizations
N = 1000             # Number of particles
dt = 0.001            # Time step
T_relax = 1000
steps_relax = int(T_relax / dt)
T = 20000              # Total simulation time
steps = int(T / dt)  # Number of steps
tau = 0.5            # Persistence time
mu = 1.0             # Mobility
ks = np.linspace(4., 10., 10)               # Spring constant
D = 100.0              # Noise strength
lb = 5.             # Trap size


# Split realizations among processes
M_local = M // size  # Number of realizations per process
if rank < M % size:
    M_local += 1



if rank == 0:
    pbar = tqdm(total=M_local * len(ks), desc="Processing", ncols=100)


res = []

for k in ks:
    # Storage for local results
    tesc_local = []

    # Realizations loop
    for m in range(M_local):
        # Initialize positions and velocities
        r = np.zeros(N)
        p = np.zeros(N)

        # Relaxation loop
        for t in range(steps_relax):
            # Add ghost particles
            r_ghost = np.concatenate(([r[0]], r, [r[-1]]))
            p_ghost = np.concatenate(([p[0]], p, [p[-1]]))

            # Calculate spring forces on r and p
            r_diff = r_ghost[2:] - 2 * r_ghost[1:-1] + r_ghost[:-2]
            p_diff = p_ghost[2:] - 2 * p_ghost[1:-1] + p_ghost[:-2]

            # Update velocities and positions using the recasted equation
            noise = np.sqrt(2 * D * dt) * np.random.randn(N) / \
                tau  # Gaussian white noise
            dp = (-p + mu * k * tau * p_diff + mu * k * r_diff) / tau
            # update r first to ensure Euler's method
            r += p * dt
            p += dp * dt + noise

        # Simulation loop
        start = 0
        for t in range(steps):
            ell = np.diff(r)
            ell_c = ell[int(N/2)-1]

            if start > 0:
                if ell_c >= lb:
                    tesc_local.append((start) * dt)
                    break
                start += 1

            if np.abs(ell_c) < 0.1 and start == 0:
                start = 1

            # Add ghost particles
            r_ghost = np.concatenate(([r[0]], r, [r[-1]]))
            p_ghost = np.concatenate(([p[0]], p, [p[-1]]))

            # Calculate spring forces on r and p
            r_diff = r_ghost[2:] - 2 * r_ghost[1:-1] + r_ghost[:-2]
            p_diff = p_ghost[2:] - 2 * p_ghost[1:-1] + p_ghost[:-2]

            # Update velocities and positions using the recasted equation
            noise = np.sqrt(2 * D * dt) * np.random.randn(N) / \
                tau  # Gaussian white noise
            dp = (-p + mu * k * tau * p_diff + mu * k * r_diff) / tau
            # update r first to ensure Euler's method
            r += p * dt
            p += dp * dt + noise

        # Store the final positions
        if rank == 0:
            pbar.update(1)

    tesc_all = comm.gather(tesc_local, root=0)

    # Combine results on the root process
    if rank == 0:
        
        ts_combined = [t for ts_process in tesc_all for t in ts_process]
        t_esc = np.array(ts_combined)

        # print(np.mean(t_esc), len(t_esc))

        tf = T * np.ones(M)
        tf[: len(t_esc)] = t_esc
        tf = t_esc
        # print(np.mean(tf), len(tf))
        print("\n", len(tf))
        res.append(np.mean(tf))

if rank == 0:
    pbar.close()
    print("\n", "Varying spring constant k:\n", res)
