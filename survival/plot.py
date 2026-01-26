import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from matplotlib import pyplot as plt



kbs = [0.005, 0.02, 0.1, 1.]  # unit: 1/h
kd = np.log(2)/16.  # unit: 1/h


start, end = 0.1, 0.92
cmap = plt.cm.plasma
colors = cmap(np.linspace(start, end, len(kbs)))



fig, ax = plt.subplots(figsize=(3, 3))

for kb in kbs:
    data = np.load(f'data/data_kb{kb}.npy')

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
    kmf.fit(T, event_observed=E, label=r"$k_b=%g\,\mathrm{h}^{-1}$" % kb, timeline=np.arange(0, 80, 0.1))

    kmf.plot_survival_function(ax=ax, color=colors[kbs.index(kb)], clip_on=False)

tt = np.arange(0, 80, 0.01)
ax.plot(tt, np.exp(-kd*tt), color="black", linestyle='--', label=r"$\mathrm{e}^{-k_d t}$")

ax.legend(frameon=False, fontsize=13)
# ax.set_yscale("log")
ax.set_xlim(0, 80)
ax.set_xticks(np.arange(0, 81, 20))

ax.set_ylim(0, 1.0)
ax.set_xlabel("Time [h]")
ax.set_ylabel("Survival probability $S(t)$")
plt.savefig("survival_prob.png", dpi=150, bbox_inches='tight')