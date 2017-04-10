# pyChirpZ

[![Build Status](https://travis-ci.org/ericmjonas/pychirpz.svg?branch=master)](https://travis-ci.org/ericmjonas/pychirpz)

Implementation of the chirp-z transform in python. Can be used to
evaluate creatively on the unit disk, or to zoom the FFT. Two
implementations, one in numba and one in C++ with eigen that is
stand-alone enough to be used in other C++ projects.

Two papers:

> Rabiner, L. R., Schafer, R. W., & Rader, C. M. (1969). The Chirp
> z-Transform Algorithm and Its Application. Bell System Technical
> Journal, 48(5),
> 1249–1292. http://doi.org/10.1002/j.1538-7305.1969.tb04268.x

> Rabiner, L. R., & Schafer, R. W. (1969). The Chirp z-Transform
> Algorithm. IEEE Transactions on Audio and Electroacoustics, 17(2),
> 86–92. http://doi.org/10.1109/TAU.1969.1162034

## Example

A basic example can be seen in [this 1D tutorial](https://github.com/ericmjonas/pychirpz/blob/master/examples/basic%20example.ipynb). 

## Installation

### Linux
Make sure you have both  [FFTW](http://www.fftw.org/) and [Boost](http://www.boost.org/)
installed. On ubuntu-based linux this can be done via
```
apt-get install  libfftw3-dev  libboost-dev
```
Note you will probably want a bleeding-edge eigency
```
pip install git+https://github.com/wouterboomsma/eigency
```

### OSX 

On OSX, make sure you have all the relevant build tools installed and install
[FFTW](http://www.fftw.org/) and [Boost](http://www.boost.org/). A recommended way is using Brew:
```
brew install boost fftw
```

Note you will probably want a bleeding-edge eigency
```
pip install git+https://github.com/wouterboomsma/eigency
```

## Benchmark 

The Chirp-Z transform lets you evaluate any evenly-spaced
set of frequencies along the unit circle (or even along an arc inside
the unit circle, but we'll ignore that right now). Imagine you 
have a 256-element-long vector, and you'd like to compute the DFT
at a more finely-spaced set of samples, but over a narrow range (the so-called
"zoomed FFT". The chirp-z transform can help. Normally we'd just
pad the FFT and then extract the region of interest in the output, 
but this can result in us doing really large FFTs. 

The speed-up can be dramatic. Below shows our input signal length (N), the number
of points in the output that we use (M), and the equivalent FFT padding size,
for a 2D fft with complex inputs. All of the numbers below are the speed
gains relative to simply padding. 

|   eq FFT points |   N |   M |   numba chirpz2d |   c++ chirpz2d32 |   c++ chirpz2d64 |
|----------------:|----:|----:|-----------------:|-----------------:|-----------------:|
|             512 |  64 |  64 |        0.548155  |         13.0742  |        10.1108   |
|             512 | 128 | 128 |        0.260806  |          3.1541  |         2.2887   |
|             512 | 256 | 256 |        0.0942346 |          1.02381 |         0.653267 |
|            1024 |  64 |  64 |        3.50962   |         82.5945  |        61.4922   |
|            1024 | 128 | 128 |        1.62635   |         16.2264  |         8.63124  |
|            1024 | 256 | 256 |        0.477032  |          4.60482 |         2.78852  |
|            2048 |  64 |  64 |       15.7544    |        345.908   |       255.567    |
|            2048 | 128 | 128 |        6.82644   |         58.3156  |        21.5589   |
|            2048 | 256 | 256 |        2.26787   |         18.6746  |         9.8219   |
|            4096 |  64 |  64 |       84.6179    |       1793.88    |      1203.37     |
|            4096 | 128 | 128 |       30.1707    |        312.881   |       223.717    |
|            4096 | 256 | 256 |        9.19736   |         87.8977  |        53.5942   |


