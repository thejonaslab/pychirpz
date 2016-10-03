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


void hello() {
    std::cout << "Hello World" << std::endl; 
    

}

ChirpZ::ChirpZ(int M, fc32_t A, fc32_t W) :
    M_(M),
    A_(A),
    W_(w)
{
    // allocate FFT 


}


ChirpZ::~ChirpZ()
{
    

}

ArrayXcf compute(ArrayXcf x) {
    
}



}
