#pragma once

#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <vector>

namespace chirpz { 

typedef std::complex<float> fc32_t;

const float PI = 3.1415926535; 
const fc32_t J(0, 1);

void hello();

class ChirpZ  {
    
public:
    ChirpZ(int M, fc32_t A, fc32_t W); 
    ~ChirpZ();
    
    ArrayXcf compute(ArrayXcf x);
    
private:
    int M_;
    fc32_t A_;
    fc32_t W_; 
}; 



}
