import numpy as np
from lifelines import KaplanMeierFitter
from matplotlib import pyplot as plt



kb = 0.02  # unit: 1/h
cell_cycles = [4, 16, 64]  # unit: h


start, end = 0.1, 0.95
cmap = plt.cm.viridis
colors = cmap(np.linspace(start, end, len(cell_cycles)))


fig, ax = plt.subplots(figsize=(3, 3))

for idx, cell_cycle in enumerate(cell_cycles):
    kd = np.log(2)/cell_cycle  # unit: 1/h

    data = np.load(f'data_kd{cell_cycle}.npy')

    thermal_time = 1000.  # time needed to reach the steady state

    timer = thermal_time
    r_times = []
    breaking = []
    arr = data
    while len(arr):
        arr = arr[arr > timer]  # Get the array after thermal_time

        if len(arr) == 0:
            break

        r_times.append(arr[0]-timer)
        breaking.append(1)
        timer += thermal_time  # Increment timer by thermal_time

    T = np.array(r_times)/kd  # Convert to hours
    E = np.array(breaking)

    print(len(T))

    kmf = KaplanMeierFitter()
    kmf.fit(T, event_observed=E,
            label=r"$k_d=%.2f\,\mathrm{h}^{-1}$" % kd, timeline=np.arange(0, 80, 0.1))

    kmf.plot_survival_function(ax=ax, color=colors[idx], clip_on=False)

    tt = np.arange(0, 80, 0.01)
    gamma = kb/kd
    eta = gamma * (np.exp(kd * tt) - 1.)
    theory = np.exp(kb*tt) * gamma * np.exp(-eta)/((1. + gamma) - np.exp(-eta))

    if idx == 2:
        ax.plot(tt, np.exp(-kd*tt), color='black', linestyle='--', label=r"$\mathrm{e}^{-k_d t}$")
    else:
        ax.plot(tt, np.exp(-kd*tt), color='black', linestyle='--')

ax.legend(frameon=False, fontsize=13)
ax.set_xlim(0, 80)
ax.set_ylim(0, 1.0)
ax.set_xlabel("Time [h]")
ax.set_xticks(np.arange(0, 81, 20))

ax.set_ylabel("Survival probability $S(t)$")

plt.savefig("survival_prob.png", dpi=150, bbox_inches='tight')
