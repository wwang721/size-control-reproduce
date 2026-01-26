import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


class Chain(object):
    def __init__(self, N: int, ratio: float = 0.1):
        self.N = N
        self.T = 0.

        self.ratio = ratio  # k_escape/k_growth

        self.Nmin = 1  # minimum chain size
        self.end_grow = False  # only the cells at the two ends can grow

        self.update_rates()

    def update_rates(self):
        growth_rate = 1.
        escape_rate = growth_rate * self.ratio

        self.k_escape = np.full(self.N-1, escape_rate)

        if self.N < 2 * self.Nmin:
            self.k_escape *= 0.
        else:
            self.k_escape[:self.Nmin-1] = np.zeros(self.Nmin-1)
            self.k_escape[self.N - self.Nmin:] = np.zeros(self.Nmin-1)

        self.k_growth = np.full(self.N, growth_rate)

        if self.end_grow and self.N > 2:
            self.k_growth[1:-1] = np.zeros(self.N - 2)

    def update(self, index):
        breaking = False
        if index < self.N-1:
            # escape
            # random pick one daughter chain
            u = np.random.randint(2)
            self.N = index + 1 if u else self.N - (index + 1)
            breaking = True
        else:
            # grow
            self.N += 1

        self.update_rates()

        return breaking

    def rfKMC(self):
        ks = np.concatenate((self.k_escape, self.k_growth))
        R = np.cumsum(ks)
        Q = np.sum(ks)

        u = 1. - np.random.uniform(0., 1.)
        index = np.searchsorted(R, u*Q)

        breaking = self.update(index)

        up = 1.-np.random.uniform(0., 1.)
        dt = np.log(1./up)/Q
        self.T += dt

        return breaking


steps = 10000000

rupture_times = []

# [0.005, 0.02, 0.1, 1.] run different kb one by one
kb = 0.005            # breaking rate
# kb = 0.02
# kb = 0.1
# kb = 1.0

kd = np.log(2)/16.  # division rate

A = Chain(500, kb/kd)

pbar = tqdm(total=steps, desc="Processing")
for _ in range(steps):
    breaking = A.rfKMC()
    if breaking:
        rupture_times.append(A.T)
    
    if _ % 1000 == 0:
        pbar.update(1000)

pbar.close()

np.save(f"data/data_kb{kb}.npy", rupture_times)


