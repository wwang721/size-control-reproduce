import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


class Chain(object):
    def __init__(self, N: int, ratio: float=0.1):
        self.N = N
        self.T = 0.
        
        self.ratio = ratio # k_escape/k_growth
        
        self.Nmin = 1 # minimum chain size
        self.end_grow = True # only the cells at the two ends can grow
        
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
        if index < self.N-1:
            # escape
            # random pick one daughter chain
            u = np.random.randint(2)
            self.N = index + 1 if u else self.N - (index + 1)
        else:
            # grow
            self.N += 1
        
        self.update_rates()

    
    def rfKMC(self):
        ks = np.concatenate((self.k_escape, self.k_growth))
        R = np.cumsum(ks)
        Q = np.sum(ks)
        
        u = 1. - np.random.uniform(0., 1.)
        index = np.searchsorted(R, u*Q)
        
        self.update(index)
        
        up = 1.-np.random.uniform(0., 1.)
        dt = np.log(1./up)/Q
        self.T += dt


N = 10000 # number of assays
steps = 100000

sizes = []
pbar = tqdm(total=N, desc="Processing")
for _ in range(N):
    A = Chain(100, 0.1)
    for _ in range(steps):
        N_tmp = A.N
        A.rfKMC()
        if A.T > 2500.:
            sizes.append(N_tmp)
            break
    pbar.update(1)
    
pbar.close()

np.save("data/2ends_k0p1", sizes)