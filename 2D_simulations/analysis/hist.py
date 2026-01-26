import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors



vas = np.arange(0.55, 1.01, 0.01)

data = np.load('fracture_rates.npy')
fractures_per_hour_connected = data[1][:len(vas)]

start, end = 0.1, 0.92
cmap = plt.cm.plasma
mycolors = cmap(np.linspace(start, end, len(vas)))


pick = [0, 5, 10, 20]
vas = [vas[i] for i in pick]
fractures_per_hour_connected = [fractures_per_hour_connected[i] for i in pick]



for idx, va in enumerate(vas):
    fig, ax = plt.subplots(figsize=(3.7, 2.5))


    data = np.load('../data/va%g.npz' % va)
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

    # ax[idx].plot(unique_states, lineage_distribution, 'x-', zorder=1)
    facecolor = colors.to_rgba(mycolors[pick[idx]], alpha=0.5)
    ax.bar(unique_states, lineage_distribution, width=1,
                align='center', color=facecolor, edgecolor='k', zorder=0)

    
    kb = fractures_per_hour_connected[idx]
    kd = np.log(2)/16.
    gamma = kb / kd
    nn = np.linspace(1, min(max(unique_states)+1, 40), 1000)
    pn = gamma/ ((1.+gamma)**nn)
    
    ax.plot(nn, pn, '--', color="k", linewidth=1.5, label=r"1D theory", zorder=2)

    ax.legend(frameon=False, fontsize=12, loc='upper right')
    ax.text(0.5, 0.75, r"$v_\mathrm{rms} =%g\,$μ$\mathrm{m/h}$" % (va*40), transform=ax.transAxes, fontsize=12, verticalalignment='top')

    ax.set_xlim(0, right=max(10, min(max(unique_states), 12)))

    if idx == 3:
        ax.set_xlim(0, 12)
        ax.set_yticks(np.arange(0, 0.9, 0.2))
    
    if idx == 0:
        ax.set_yticks(np.arange(0, 0.6, 0.1))
        
    ax.set_xticks(np.arange(2, 13, 2))
    ax.set_ylabel("$p_n$")

    ax.set_xlabel("$n$")


    # plt.subplots_adjust(wspace=0.3)  # Increase spacing
    # plt.show()
    plt.savefig(f"hist{idx}.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
