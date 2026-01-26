import numpy as np
import matplotlib.pyplot as plt


def matrix_method(k1, k2, Nmin, order):
    # Left-hand side
    A = np.zeros((order, order))

    for i in range(order):
        A[i, i] = -k1 * np.maximum(0, Nmin+i - 2 * Nmin + 1) - k2 * (Nmin+i)

    for i in range(1, order):
        A[i, i-1] = k2 * (Nmin+i - 1)

    order2 = order - Nmin  # order for small upper right matrix
    B = np.triu(k1*np.ones((order2, order2)))

    A[:order2, -order2:] += B

    C = np.ones(order)  # normalization

    K = np.zeros((order+1, order))
    K[:-1] += A
    K[-1] += C

    # Right-hand side
    D = np.zeros(order+1)
    D[-1] = 1.

    # Compute the solution
    res = np.linalg.inv(K.T @ K) @ K.T @ D.T

    return res



#-------------------------------------------------------------
# Define the range you want for the new colormap (0.1 to 0.95)
start, end = 0.1, 0.95
cmap = plt.cm.viridis
colors = cmap(np.linspace(start, end, 4))

fig, ax2 = plt.subplots(1, 1, figsize=(3.5, 3), dpi=150)

Nmin = 200  # minimum cluster size

size = 5000  # matrix size
nn = np.arange(Nmin, Nmin+size)


gamma = 0.01 # the ratio of k_esc/k_grow 
k2 = 1.
k1 = k2 * gamma
res = matrix_method(k1, k2, Nmin, size)
ax2.plot(nn, res, '--', color=colors[0], label=r"$\gamma=0.01$", linewidth=2)

sizes = np.load('min200_k0p01.npy')
interval = 20
bins = np.arange(200, 1500, interval)  # From 10 to 30 with interval 5
hist, bin_edges = np.histogram(sizes, bins=bins, density=True)
bin_center = bin_edges[:-1] + (interval/2.)
ax2.plot(bin_center[bin_center <= 800], hist[bin_center <= 800], "x", color=colors[0], clip_on=False, markeredgewidth=1.2)


gamma = 0.1 # the ratio of k_esc/k_grow 
k2 = 1.
k1 = k2 * gamma
res = matrix_method(k1, k2, Nmin, size)

ax2.plot(nn, res, color=colors[1], label=r"$\gamma=0.1$", linewidth=2)

sizes = np.load('min200_k0p1.npy')
interval = 20
bins = np.arange(200, 800, interval)  # From 10 to 30 with interval 5
hist, bin_edges = np.histogram(sizes, bins=bins, density=True)
bin_center = bin_edges[:-1] + (interval/2.)
ax2.plot(bin_center[hist!=0], hist[hist!=0], "x", color=colors[1], clip_on=False, markeredgewidth=1.2)


ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax2.set_xticks(np.linspace(300, 700, 3), minor=True)
ax2.set_xlim(180, 800)
ax2.set_ylim(bottom=0)
ax2.set_xlabel(r"$n$", fontsize=16)
ax2.set_ylabel(r"$p_n$", fontsize=16)
ax2.legend(frameon=False)

plt.savefig('min_N.png', bbox_inches='tight')
