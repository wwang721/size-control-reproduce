import numpy as np
from tqdm import tqdm
from mpi4py import MPI

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Parameters
M = 48 * 100        # Total number of realizations
N = 1000        # Number of particles
dt = 0.01       # Time step

# Need a smaller dt for smaller tau
# dt = 0.0001        # smaller time step

T = 1000       # Total simulation time
steps = int(T / dt)     # Number of steps
taus = 10**np.linspace(-2, 2, 20)     # Persistence time
mu = 1.0        # Mobility
k = 1.0     # Spring constant
D = 1.0     # Noise strength

# Split realizations among processes
M_local = M // size  # Number of realizations per process
if rank < M % size:
    M_local += 1


if rank == 0:
    pbar = tqdm(total=M_local*len(taus), desc="Processing")

index = 0
for tau in taus:
    # Storage for local results
    rs_local = []

    # Process realizations
    for m in range(M_local):
        # Initialize positions and velocities
        r = np.zeros(N)
        p = np.zeros(N)

        # Simulation loop
        for t in range(steps):
            # Add ghost particles
            r_ghost = np.concatenate(([r[0]], r, [r[-1]]))
            p_ghost = np.concatenate(([p[0]], p, [p[-1]]))

            # Calculate spring forces on r and p
            r_diff = r_ghost[2:] - 2 * r_ghost[1:-1] + r_ghost[:-2]
            p_diff = p_ghost[2:] - 2 * p_ghost[1:-1] + p_ghost[:-2]

            # Update velocities and positions using the recasted equation
            noise = np.sqrt(2 * D * dt) * np.random.randn(N) / tau  # Gaussian white noise
            dp = (-p + mu * k * tau * p_diff + mu * k * r_diff) / tau
            # update r first to ensure Euler's method
            r += p * dt
            p += dp * dt + noise

        # Store the final positions
        rs_local.append(r)
        if rank == 0:
            pbar.update(1)

    # Gather results from all processes
    rs_all = comm.gather(rs_local, root=0)

    # Combine results on the root process
    if rank == 0:
        rs_combined = [r for rs_process in rs_all for r in rs_process]
        rs_combined = np.array(rs_combined)

        # save results
        np.save("data/rs_%d" % index, rs_combined)
        index += 1

if rank == 0:
    pbar.close()
