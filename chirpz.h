#pragma once

#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <vector>

namespace chirpz { 

using namespace Eigen; 
typedef std::complex<float> fc32_t;

const float PI = 3.1415926535; 
const fc32_t J(0, 1);

void hello(fc32_t x);

class ChirpZ  {
    
public:
    ChirpZ(int N, int M, fc32_t A, fc32_t W); 
    ~ChirpZ(void);
    
    ArrayXcf compute(const ArrayXcf & x);
    
private:
    int N_; 
    int M_;
    fc32_t A_;
    fc32_t W_;

    int L_;
    ArrayXcf yn_scale_; 

    fftwf_complex* fftw_in_;
    fftwf_complex* fftw_out_;
    
    fftwf_plan fftw_forward_plan_;
    fftwf_plan fftw_reverse_plan_;

    ArrayXcf Vr_;
    ArrayXcf g_scale_;

    // convenience functions, warning, make copies
    ArrayXcf fft(const ArrayXcf & in);
    ArrayXcf ifft(ArrayXcf in);
        

}; 


class ChirpZ2d  {
    
public:
    ChirpZ2d(int N, int M, fc32_t A, fc32_t W); 
    ~ChirpZ2d(void);
    
    MatrixXcf compute(const MatrixXcf & x);
    
private:
    int N_; 
    int M_;
    fc32_t A_;
    fc32_t W_;

    int L_;
    MatrixXcf yn_scale_; 

    fftwf_complex* fftw_in_;
    fftwf_complex* fftw_out_;
    
    fftwf_plan fftw_forward_plan_;
    fftwf_plan fftw_reverse_plan_;

    MatrixXcf Vr_;
    MatrixXcf g_scale_;

    // convenience functions, warning, make copies
    MatrixXcf fft(const MatrixXcf & in);
    MatrixXcf ifft(const MatrixXcf & in);
        

}; 



}
