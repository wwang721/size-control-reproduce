import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


vas = np.arange(0.55, 1.51, 0.01)


start, end = 0.1, 0.92
cmap = plt.cm.plasma
colors = cmap(np.linspace(start, end, len(vas[:58])))



means = []
stds = []
for va in vas:
    data = np.load('../data/va%g.npz' % va)
    tot_num_breaking = data['breaking']
    tot_num_connection_hours = data['connections']

    fracture_per_hour_connected = tot_num_breaking / tot_num_connection_hours

    means.append(np.mean(fracture_per_hour_connected))
    stds.append(np.std(fracture_per_hour_connected))

fig, ax = plt.subplots(figsize=(3.6, 3))

print(vas.tolist())
print(np.array(means).tolist())
print(np.array(stds).tolist())
np.save('fracture_rates.npy', np.array([vas, means, stds]))


kb_theory = np.linspace(0.6, 1.1, 100)
kd = np.log(2)/16.
gamma_theory = kb_theory / kd
mean_theory = (1. + gamma_theory) / gamma_theory
var_theory = mean_theory / gamma_theory
CV_theory = 1./np.sqrt(1.+gamma_theory)


vv = np.linspace(0.2, 1.5, 100)
muk = 4 # 1/4 [hr^-1]
deltalb = 6.5 # [um]
tau = 0.5 # [hr]
Ds = (vv * 40)**2 * tau / 2.
eta = 1. + 2. * muk * tau
tau0 = np.pi * eta / muk
kb_1d = (1./tau0) * np.exp(-(muk * (deltalb**2) * np.sqrt(1. + 4 * muk * tau))/(2 * Ds))


for idx in range(58):
    ax.errorbar((vas * 40)[idx], means[idx], yerr=stds[idx], fmt='o', color=colors[idx], markerfacecolor='None', markeredgewidth=1.5, capsize=4, elinewidth=1.5, clip_on=False, zorder=0)

ax.set_xlabel("$v_\mathrm{rms}$ [μ$\mathrm{m/h}$]")
ax.set_ylabel("Fractures per hour connected")

ax.text(23, 1.25, f"SD from {len(fracture_per_hour_connected)} simulations", fontsize=12)



# ax.set_yscale("log")
ax.set_xlim(0.54 * 40, 45)
ax.set_ylim(0, 1.4)
# ax.set_yticks(np.arange(0, 1.6, 0.5))
ax.set_xticks(np.arange(25, 46, 5))

ax.yaxis.set_minor_locator(AutoMinorLocator(2))


plt.savefig("fracture_rate.png", dpi=150, bbox_inches='tight')
