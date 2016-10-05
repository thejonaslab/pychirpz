import numpy as np
import pychirpz


def test_fft_eq():
    
    N = 128

    for M in [128, 32, 10]:
        for start_idx in [0, 4, 8]:

            x = np.random.normal(0, 1, N)

            x_fft = np.fft.fft(x)

            w_delta = 2.0*np.pi/N

            start = start_idx * w_delta
            x_chirpz = pychirpz.zoom_fft(x, start, w_delta, M)

            np.testing.assert_allclose(np.roll(x_fft, -start_idx)[:M],
                                       x_chirpz[:M], verbose=True)


def test_fft2d_eq():
    
    N = 256
    for M in [20, 128, N]:
        for start_idx in [0, 5, -5]:

            x = np.random.normal(0, 1, (N, N))


            x_fft = np.fft.fft2(x)

            w_delta = 2.0*np.pi/N

            start = start_idx * w_delta
            x_chirpz = pychirpz.zoom_fft2(x, start, w_delta, M)

            x_fft = np.roll(x_fft, -start_idx, axis=0)
            x_fft = np.roll(x_fft, -start_idx, axis=1)

            np.testing.assert_allclose(x_fft[:M, :M],
                                       x_chirpz[:M, :M], verbose=True)
