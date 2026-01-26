import numpy as np
import scipy.special
import matplotlib.pyplot as plt


plt.figure(figsize=(4,3), dpi=200)
ax1 = plt.subplot(1,1,1)

gamma = 0.1

nn = np.arange(1, 101)
pn_all = gamma/((1+gamma)**nn)

ax1.plot(nn, pn_all, color="C0", label=r"All", linewidth=2, zorder=1)

sizes = np.load('data/all_k0p1.npy')
unique_integers, counts = np.unique(sizes, return_counts=True)
total_count = np.sum(counts)
probabilities = counts / total_count
ax1.plot(unique_integers[:20], probabilities[:20], "x", clip_on=False, zorder=1)

#============================================

nn = np.linspace(2, 25, 100)
nn = np.insert(nn, 0, 1)
theory = nn * ((2./gamma)**(nn-1))/scipy.special.poch((2.+gamma)/gamma, nn)
ax1.plot(nn, theory, '-.', color="C4", label="One-end", linewidth=2, zorder=2)

sizes = np.load('data/1end_k0p1.npy')
unique_integers, counts = np.unique(sizes, return_counts=True)
total_count = np.sum(counts)
probabilities = counts / total_count
ax1.plot(unique_integers[:20], probabilities[:20], "x", color="C4", clip_on=False, zorder=1)

#============================================

theory *= (2.+gamma)/(2 * (1+gamma))
theory[0] = gamma/(1.+gamma)
plt.plot(nn, theory, '--', color="C3", label="Two-end", linewidth=2, zorder=3)

sizes = np.load('data/2ends_k0p1.npy')
unique_integers, counts = np.unique(sizes, return_counts=True)
total_count = np.sum(counts)
probabilities = counts / total_count
ax1.plot(unique_integers[:20], probabilities[:20], "x", color="C3", clip_on=False, zorder=3)


ax1.set_xlim(1, 20)
ax1.set_xlabel(r"$n$")
ax1.legend(frameon=False)
ax1.set_ylabel(r"$p_n$")

plt.savefig("distributions.png", bbox_inches='tight')