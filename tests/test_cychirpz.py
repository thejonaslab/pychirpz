import unittest
import numpy as np
import sys
from chirpz import cychirpz
import pyfftw

# import pyximport; 
# pyximport.install()
# import cychirpz
# import pychirpz
# import pyfftw


#@unittest.skip("temporarily disabled")
def test_fft_eq():
    
    N = 128

    for M in [10, 20, 128]:
        for start_idx in [0, 5]:

            x = np.random.normal(0, 1, N)

            x_fft = np.fft.fft(x)
            x_fft = np.roll(x_fft, -start_idx)[:M]
            w_delta = 2.0*np.pi/N

            start = start_idx * w_delta
            x_chirpz = cychirpz.zoom_fft(x, start, w_delta, M)

            # for i in range(M):
            #     print i, np.abs(x_fft[i] - x_chirpz[i]), np.abs(x_fft[i])

            np.testing.assert_allclose(x_fft[:M],
                                       x_chirpz[:M], verbose=True, 
                                       atol=1e-5, rtol=1e-5)


#@unittest.skip("temporarily disabled")
def test_fft2d_eq():
    
    N = 256
    for M in [10, 20, 256]:
        for start_idx in [0, 5]:

            x = np.random.normal(0, 1, (N, N)) # .astype(np.float32)

            x_fft = pyfftw.interfaces.numpy_fft.fft2(x).T

            w_delta = 2.0*np.pi/N

            start = start_idx * w_delta

            x_chirpz = cychirpz.zoom_fft2(x, start, w_delta, M)

            print(x_chirpz.shape)

            x_fft = np.roll(x_fft, -start_idx, axis=0)
            x_fft = np.roll(x_fft, -start_idx, axis=1)

            np.testing.assert_allclose(x_fft[:M, :M],
                                       x_chirpz[:M, :M], verbose=True, 
                                       rtol=1e-7) # FIXME WOW THIS IS LOW
