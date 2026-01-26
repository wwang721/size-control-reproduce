import numpy as np
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter



T_tot = 1000.  # Total time of simulation


vas = [0.55, 0.6, 0.8, 1]

start, end = 0.1, 0.8
cmap = plt.cm.plasma
colors = cmap(np.linspace(start, end, len(vas)))

fig, ax = plt.subplots(figsize=(3.6, 3))

for va in vas:
    data = np.load('../data/va%g.npz' % va)

    rupture_times = data['rupture']

    # Find indices where rupture_times == 0
    zero_idx = np.where(rupture_times == 0.)[0]

    # Add end of array as final boundary
    zero_idx = np.append(zero_idx, len(rupture_times))

    # Extract segments between the 0s
    segments = [rupture_times[zero_idx[i]+1: zero_idx[i+1]]
                for i in range(len(zero_idx)-1)]

    Ts = []
    Es = []
    for segment in segments:
        thermal_time = 0.2 * T_tot  # time needed to reach the steady state

        timer = thermal_time
        r_times = []
        breaking = []
        arr = segment
        while timer < T_tot:
            arr = arr[arr > timer]  # Get the array after thermal_time
            if len(arr) != 0:
                r_times.append(arr[0]-timer)
                breaking.append(1)
                timer += thermal_time  # Increment timer by thermal_time
            else:
                r_times.append(T_tot - timer)
                breaking.append(0)
                break

        T = np.array(r_times)
        E = np.array(breaking)

        Ts.append(T)
        Es.append(E)

    T = np.concatenate(Ts)
    E = np.concatenate(Es)

    # print(len(E))

    kmf = KaplanMeierFitter()
    kmf.fit(T, event_observed=E,
            label=r"$v_\mathrm{rms}=%g\,$μ$\mathrm{m/h}$" % (va*40), timeline=np.arange(0, 100.1, 0.1))

    idx = vas.index(va)
    kmf.plot_survival_function(ax=ax, color=colors[idx], clip_on=False)

kd = np.log(2)/16.  # unit: 1/h
tt = np.arange(0, 100, 0.01)
ax.plot(tt, np.exp(-kd*tt), color='black', linestyle='--', label=r"$\mathrm{e}^{-k_d t}$")

ax.legend(frameon=False, fontsize=10)
# ax.set_yscale("log")
ax.set_xlim(0, 100)
ax.set_ylim(0, 1.0)
ax.set_xlabel("Time [h]")
ax.set_ylabel("Survival probability $S(t)$")


plt.savefig("survival.png", dpi=150, bbox_inches='tight')
