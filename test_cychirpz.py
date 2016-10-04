import numpy as np
import sys
import pyximport; 
pyximport.install()
import cychirpz



N = 128
M = 10
np.random.seed(0)

x = np.random.normal(0, 1, N)

x_fft = np.fft.fft(x)

w_delta = 2.0*np.pi/N
start_idx = 0
start = start_idx * w_delta

A = np.exp(1j * start)
W = np.exp(-1j * w_delta)

cz = cychirpz.PyChirpZ(N, M, A, W)
x = x.astype(np.complex64)

a = cz.compute(x)
print a[:5]
print x_fft[:5]
