import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator



vas = np.arange(0.55, 1.51, 0.01)

data = np.load('fracture_rates.npy')
kbs= data[1][:len(vas)]
stds = data[2][:len(vas)]



vas = vas[::2]
kbs = kbs[::2]
stds = stds[::2]

start, end = 0.1, 0.92
cmap = plt.cm.plasma
colors = cmap(np.linspace(start, end, len(vas)))


means = []
vars = []
gammas_fitted = []

for va in vas:
    data = np.load(r'../data/va%g.npz' % va)
    tt = data['t']
    N_trajectory = data['N']

    start_idx = int(0.2 * len(tt))
    # Skip the first 20% of the data
    N_trajectory = N_trajectory[:, start_idx:]

    # Step 1: Compute time intervals
    dt = np.diff(tt)[-1]  # Assuming uniform time intervals
    sizes = np.concatenate(N_trajectory)  # Cluster sizes at each time step

    # Step 2: Compute total time spent in each state using NumPy
    unique_states, inverse_indices = np.unique(sizes[:-1], return_inverse=True)
    time_in_state = np.zeros_like(unique_states, dtype=float)

    # Accumulate total time for each state
    np.add.at(time_in_state, inverse_indices, dt)

    # Step 4: Normalize to get the distribution
    lineage_distribution = time_in_state / np.sum(time_in_state)

    avg = np.sum(unique_states * lineage_distribution)
    variance = np.sum(((unique_states - avg) ** 2) * lineage_distribution)


    means.append(avg)
    vars.append(variance)




fig_size = (3.5, 3)

fig, ax = plt.subplots(figsize=fig_size)

kb_theory = 10**np.linspace(-2, 3, 100)
kd = np.log(2)/16.
gamma_theory = kb_theory / kd
mean_theory = (1. + gamma_theory) / gamma_theory
var_theory = mean_theory / gamma_theory
CV_theory = 1./np.sqrt(1.+gamma_theory)

kbs = np.array(kbs)
stds = np.array(stds)/kd

for idx in range(len(vas)):
    if idx == 0:
        ax.errorbar((kbs/kd)[idx], means[idx], xerr=stds[idx], color=colors[idx], fmt='o', label="Simulations", markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)
    else:
        ax.errorbar((kbs/kd)[idx], means[idx], xerr=stds[idx], color=colors[idx], fmt='o', markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)

ax.plot(kb_theory/kd, mean_theory, 'k--', label='1D theory', zorder=2)
ax.set_xscale('log')
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel(r"$\langle n\rangle$")
ax.set_xlim(1, 70)
ax.set_ylim(1, 2.4)
ax.set_yticks(np.arange(1.2, 2.5, 0.4))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.legend(frameon=False)
ax.tick_params(axis='x', pad=5)
plt.savefig("mean.png", dpi=150, bbox_inches='tight')
plt.close()


fig, ax = plt.subplots(figsize=fig_size)

for idx in range(len(vas)):
    if idx == 0:
        ax.errorbar((kbs/kd)[idx], vars[idx], xerr=stds[idx], color=colors[idx], fmt='o', label="Simulations", markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)
    else:
        ax.errorbar((kbs/kd)[idx], vars[idx], xerr=stds[idx], color=colors[idx], fmt='o', markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)

ax.plot(kb_theory/kd, var_theory, 'k--', label='1D theory', zorder=2)
# ax.plot(kbs_fitted, vars, '^', color="C3", label="Fitted $k_b$", markerfacecolor='None', markeredgewidth=1.5, clip_on=False, zorder=3)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel(r"$\langle n^2\rangle - \langle n\rangle^2$")
ax.set_xlim(1, 70)
ax.set_ylim(2e-2, 5)
ax.legend(frameon=False)
ax.tick_params(axis='x', pad=5)
plt.savefig("variance.png", dpi=150, bbox_inches='tight')
plt.close()


fig, ax = plt.subplots(figsize=fig_size)


ax.set_xscale('log')

for idx in range(len(vas)):
    if idx == 0:
        ax.errorbar((kbs/kd)[idx], (np.sqrt(vars)/np.array(means))[idx], xerr=stds[idx], color=colors[idx], fmt='o', label="Simulations", markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)
    else:
        ax.errorbar((kbs/kd)[idx], (np.sqrt(vars)/np.array(means))[idx], xerr=stds[idx], color=colors[idx], fmt='o', markerfacecolor='None', markeredgewidth=1.5, capsize=4, zorder=1)

ax.plot(kb_theory/kd, CV_theory, 'k--', label='1D theory', zorder=2)
# ax.plot(kbs_fitted, np.sqrt(vars)/np.array(means), '^', label="Fitted $k_b$", color="C3", markerfacecolor='None', markeredgewidth=1.5, clip_on=False, zorder=3)
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel(r"$\mathrm{CV}=\sigma_n/\langle n\rangle$")
ax.set_xlim(1, 70)
ax.set_ylim(0.12, 1)
ax.legend(frameon=False)
ax.tick_params(axis='x', pad=5)
plt.savefig("CV.png", dpi=150, bbox_inches='tight')

