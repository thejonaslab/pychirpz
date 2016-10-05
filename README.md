# pyChirpZ

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

You have a function defined with finite support and you wish to
evaluate the the DFT on it. 


## Benchmarks


## code organization


## To Do 
- [ ] make real setup.py, cython build
- [ ] Find a way to make cython build optional 
- [ ] Create templatized versions of functions to allow float/double operations
- [ ] Compute the complex exponentiations in a more numerically-stable way 
- [ ] We should be using block operations for eigen instead of iterating 
- [ ] More checks for types passed into eigen functions
- [ ] Fix the insane row/col major errors
