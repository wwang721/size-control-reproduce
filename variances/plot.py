import numpy as np
import matplotlib.pyplot as plt

taus = 10**np.linspace(-2, 2, 20)
N = 1000 # number of particles

res = []

index = 0
for tau in taus:
    rs = np.load("data/rs_%d.npy" % index)
    ells = np.diff(rs, axis=1)
    
    means = []
    variances = []
    for it in range(N-1):
        ll = ells[:, it]
        means.append(np.mean(ll))
        variances.append(np.var(ll))

    res.append(variances[int(N/2)-1])
    index += 1
    
print(res)

simulations = np.array(res)

# params:
D = 1.    # diffusion
mu = 1.   # motility
k = 1.    # spring constant

tau = 10**np.linspace(-2, 2, 100)

theory_DFT = D/(mu*k*np.sqrt(1.+4*mu*k*tau))

theory_asymptot = 1. / (2. * np.sqrt(mu * k * tau))
theory_2particle = 1. / (1. + 2. * mu * k * tau)


plt.figure(figsize=(3.8, 3), dpi=150)
ax = plt.subplot(1,1,1)

# plot theories
ax.axhline(1., color="#A9A9A9", linestyle="dashed", linewidth=2., zorder=0)
ax.plot(tau, theory_DFT, color="C0", label=r"$1∕\sqrt{1+4\mu k\tau}$", linewidth=2., zorder=1)
ax.plot(tau, theory_2particle, "--", color="C3", label=r"$1∕(1+2\mu k\tau)$", linewidth=2., zorder=2)

# plot simulations
tau = 10**np.linspace(-2, 2, 20)
ax.plot(tau, simulations, 'o', color="C0", label=r"Simulations", markeredgewidth=1.5, markerfacecolor="None", clip_on=False, zorder=4)

ax.set_xscale('log')
ax.set_ylim(0.001, 1.15)
ax.set_xlim(0.01, 100)

plt.xlabel(r"$\mu k\tau$")
plt.ylabel(r"$\langle\ell_n^2\rangle/\langle\ell^2\rangle_\mathrm{th}$")
plt.legend(frameon=False, bbox_to_anchor=(0.50, 0.38), fontsize=10)
plt.savefig("variances.png", bbox_inches='tight')
