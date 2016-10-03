import numpy as np
import pychirpz


def test_fft_eq():
    
    N = 128
    M = 10
    x = np.random.normal(0, 1, N)

    x_fft = np.fft.fft(x)
    
    w_delta = 2.0*np.pi/N
    start_idx = 40
    start = start_idx * w_delta
    x_chirpz = pychirpz.zoom_fft(x, start, w_delta, M)

    np.testing.assert_allclose(np.roll(x_fft, -start_idx)[:M],
                               x_chirpz[:M], verbose=True)
