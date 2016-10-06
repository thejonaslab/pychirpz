#include <fftw3.h>
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>
#include <cmath>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <unsupported/Eigen/FFT>
#include <fstream>
#include "chirpz.h"


using namespace Eigen; 
namespace chirpz {


void hello(fc32_t c) {
    std::cout << "Hello World : " << c << std::endl; 
    

}

int nextpo2(int x) {
    return int(std::ceil(std::log2(x)));
    
}

template class ChirpZ<c32_t>;
template class ChirpZ<c64_t>;
template class ChirpZ2d<c32_t>;
template class ChirpZ2d<c64_t>;

}

