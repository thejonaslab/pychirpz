import numpy as np
import time
import chirpz
import chirpz.cychirpz
import pandas as pd
import tabulate
import cPickle as pickle
from ruffus import * 

ITERS = 20


@files(None, "benchmark.pickle")
def run(infile, outfile):

    rows = []

    for FFT_N in [512, 1024, 2048, 4096]:
        for N in [64, 128, 256, 512]:
            for M in [32, 64, 128, 256]:
                print FFT_N, N, M

                x2d = np.zeros((N, N))
                width = 8
                center = N // 2
                x2d[center-width:center+width, center-width:center+width] = 1.0

                fft_df = 2*np.pi / FFT_N 
                phi_0 = fft_df / (2*np.pi)
                W = np.exp(-1.0j * 2.0*np.pi * phi_0)
                A = np.exp(-1j * np.pi/16)

                cyCZ32 = chirpz.cychirpz.PyChirpZ2d32(N, M, A, W)
                cyCZ64 = chirpz.cychirpz.PyChirpZ2d64(N, M, A, W)

                for i in range(ITERS):
                    res = {'iter' : i, 
                           'df' : fft_df, 
                           'fft N' : FFT_N, 
                           'N' : N, 
                           'M' : M}
                    fft2_t1 = time.time()
                    out2dfft = np.fft.fft2(x2d, (FFT_N, FFT_N))
                    fft2_t2 = time.time()
                    res['fft2'] = fft2_t2 - fft2_t1

                    t1 = time.time()
                    out2d = chirpz.chirpz2d(x2d, M, A, W)
                    t2 = time.time()
                    res['numba chirpz2d'] = t2-t1

                    t1 = time.time()
                    out2d = cyCZ32.compute(x2d.astype(np.complex64))
                    t2 = time.time()
                    res['c++ chirpz2d32'] = t2-t1

                    t1 = time.time()
                    out2d = cyCZ64.compute(x2d.astype(np.complex128))
                    t2 = time.time()
                    res['c++ chirpz2d64'] = t2-t1

                    if i > 0: # skip first for numba warm-up

                        rows.append(res)

    benchdf = pd.DataFrame(rows)
    pickle.dump(benchdf, open(outfile, 'w'))


if __name__ == "__main__":
    pipeline_run([run])
