import numpy as np
import matplotlib.pyplot as plt

# Change trape size lb

lbs = np.linspace(5., 10., 100)            # Trap size
tau = 0.5            # Persistence time
mu = 1.0             # Mobility
k = 4.              # Spring constant
D = 100.              # Noise strength

simulations0 = np.load('change_lb/res.npy')

lth = np.sqrt(2*D/(mu*k*np.sqrt(1.+4*mu*k*tau)))
prefactor = 2 * np.pi * (1.+2*mu*k*tau) / (2. * mu * k)
theories0 = prefactor * np.exp((lbs/lth)**2)

#====================================
# Change spring constant k

lb = 5.           # Trap size
tau = 0.5            # Persistence time
mu = 1.0             # Mobility
ks = np.linspace(4, 10, 100)              # Spring constant
D = 100.              # Noise strength

simulations1 = np.load('change_k/res.npy')

lth = np.sqrt(2*D/(mu*ks*np.sqrt(1.+4*mu*ks*tau)))
prefactor = 2 * np.pi * (1+2*mu*ks*tau) / (2. * mu * ks)
theories1 = prefactor * np.exp((lb/lth)**2)

#====================================
# Change persistence time tau

lb = 5.           # Trap size
taus = np.linspace(0.5, 4, 100)            # Persistence time
mu = 1.0             # Mobility
k = 4.              # Spring constant
D = 100.              # Noise strength

simulations2 = np.load('change_tau/res.npy')

lth = np.sqrt(2*D/(mu*k*np.sqrt(1.+4*mu*k*taus)))
prefactor = 2 * np.pi * (1+2*mu*k*taus) / (2. * mu * k)
theories2 = prefactor * np.exp((lb/lth)**2)

#====================================
# Change diffusion coefficient D

lb = 5.           # Trap size
tau = 0.5        # Persistence time
mu = 1.0             # Mobility
k = 4.              # Spring constant
Ds = np.linspace(25, 100, 100)              # Noise strength

simulations3 = np.load('change_D/res.npy')

lth = np.sqrt(2*Ds/(mu*k*np.sqrt(1.+4*mu*k*tau)))
prefactor = 2 * np.pi * (1+2*mu*k*tau) / (2. * mu * k)
theories3 = prefactor * np.exp((lb/lth)**2)

#====================================
# Plotting


fig = plt.figure(figsize=(4, 3), dpi=150)

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

ax1.plot(lbs, theories0, color="C0", linewidth=1.8)
ax1.plot(lbs[::11], simulations0, "o", color="C0", markeredgewidth=1.5, markerfacecolor="None", clip_on=False)
ax1.set_ylabel(r"$\tau_\mathrm{esc}\:[\mathrm{h}]$", fontsize=15, labelpad=2)
ax1.set_xlabel(r"$\delta\ell_{\!b}\:[$μ$\mathrm{m}]$", fontsize=15, labelpad=3)
ax1.tick_params(axis='both', labelsize=12)
ax1.set_xlim(5, 10)


ax2.plot(ks, theories1, color="C1", linewidth=1.8)
ax2.plot(ks[::11], simulations1, "o", color="C1", markeredgewidth=1.5, markerfacecolor="None", clip_on=False)
ax2.set_xlabel(r"$\mu k\:[\mathrm{h}^{-1}]$", fontsize=15, labelpad=2)
ax2.tick_params(axis='both', labelsize=12)
ax2.set_xlim(4, 10)
ax2.set_xticks([5, 7.5, 10])


ax3.plot(taus, theories2, color="C2", linewidth=1.8)
ax3.plot(taus[::11], simulations2, "o", color="C2", markeredgewidth=1.5, markerfacecolor="None", clip_on=False)
ax3.set_ylabel(r"$\tau_\mathrm{esc}\:[\mathrm{h}]$", fontsize=15, labelpad=2)
ax3.set_xlabel(r"$\tau\:[\mathrm{h}]$", fontsize=15, labelpad=4)
ax3.tick_params(axis='both', labelsize=12)
ax3.set_xlim(0.5, 4)
ax3.set_xticks([1, 2, 3, 4])


ax4.plot(Ds, theories3, color="C3", linewidth=1.8)
ax4.plot(Ds[::11], simulations3, "o", color="C3", markeredgewidth=1.5, markerfacecolor="None", clip_on=False)
ax4.set_xlabel(r"$D\:[$μ$\mathrm{m^2\!∕h]}$", fontsize=15, labelpad=2)
ax4.tick_params(axis='both', labelsize=12)
ax4.set_xlim(25, 100)
ax4.set_xticks([25, 50, 75, 100])

ax1.set_yscale("log")
ax2.set_yscale("log")
ax3.set_yscale("log")
ax4.set_yscale("log")


plt.subplots_adjust(wspace=0.3, hspace=0.6)
plt.savefig('MFPT.png', bbox_inches='tight')