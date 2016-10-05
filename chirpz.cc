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

ChirpZ::ChirpZ(int N, int M, fc32_t A, fc32_t W) :
    N_(N), M_(M), A_(A), W_(W),
    L_(std::pow(2, nextpo2(N + M -1)))
{

    // allocate FFT 
    fftw_in_ = (fftwf_complex*)fftwf_alloc_complex(L_);
    fftw_out_ = (fftwf_complex*)fftwf_alloc_complex(L_);

    fftw_forward_plan_ = fftwf_plan_dft_1d(L_, fftw_in_,
                                           fftw_out_,
                                           FFTW_FORWARD, FFTW_MEASURE);
    
    fftw_reverse_plan_ = fftwf_plan_dft_1d(L_, fftw_in_,
                                           fftw_out_,
                                           FFTW_BACKWARD, FFTW_MEASURE);
    
    yn_scale_ = ArrayXcf::Zero(L_);
    for(int n = 0; n < N_; n++) {
        yn_scale_(n) = std::pow(A_, -n) * std::pow(W_, (n*n)/2.0f); 
    }

    ArrayXcf vn = ArrayXcf::Zero(L_);
    for(int n = 0; n < M_; ++n) {
        vn(n) = std::pow(W, (-n*n)/2.0f);
    }
    
    for(int n = L_-N_+1; n < L_; ++n) {
        vn(n) = std::pow(W_, -((L_-n)*(L_-n))/2.0f); 
    }

    Vr_ = fft(vn);
            
    g_scale_ = ArrayXcf::Zero(M_);
    for(int k = 0; k < M_; ++k) {
        g_scale_(k) = std::pow(W_, (k*k/2.0f)); 
    }
        

    
}


ChirpZ::~ChirpZ(void )
{
    fftwf_free(fftw_in_);
    fftwf_free(fftw_out_); 
    

}

ArrayXcf ChirpZ::compute(const ArrayXcf & x) {


    ArrayXcf yn = ArrayXcf::Zero(L_);
    
    for(int i = 0; i < x.size(); ++i) { 
        yn(i) = x(i) * yn_scale_(i); 
    }

    
    ArrayXcf Yr = fft(yn);

    ArrayXcf Gr = Yr * Vr_;


    ArrayXcf gk = ifft(Gr);
    
    ArrayXcf Xk = ArrayXcf::Zero(g_scale_.size());
    for (int i = 0; i < g_scale_.size(); ++i) {
        Xk(i) = g_scale_(i) * gk(i);
    }
    ArrayXcf out =  Xk / (2*N_);

    return out; 

}

ArrayXcf ChirpZ::fft(const ArrayXcf & in) {
    Map<ArrayXcf> in_f((fc32_t *)fftw_in_, L_);

    in_f = in;

    fftwf_execute(fftw_forward_plan_);
    
    Map<ArrayXcf> out_field((fc32_t*) fftw_out_, L_) ;
    ArrayXcf out = out_field;
    return out; 
}

ArrayXcf ChirpZ::ifft(ArrayXcf in) {

    Map<ArrayXcf> in_f((fc32_t *)fftw_in_, L_);
    in_f = in;
    fftwf_execute(fftw_reverse_plan_);
    
    Map<ArrayXcf> out_field((fc32_t*) fftw_out_, L_) ;
    ArrayXcf out = out_field;
    return out;
    

}



ChirpZ2d::ChirpZ2d(int N, int M, fc32_t A, fc32_t W) :
    N_(N), M_(M), A_(A), W_(W),
    L_(std::pow(2, nextpo2(N + M -1)))
{

    // allocate FFT 
    fftw_in_ = (fftwf_complex*)fftwf_alloc_complex(L_ * L_);
    fftw_out_ = (fftwf_complex*)fftwf_alloc_complex(L_ * L_);

    fftw_forward_plan_ = fftwf_plan_dft_2d(L_, L_, fftw_in_,
                                           fftw_out_,
                                           FFTW_FORWARD, FFTW_MEASURE);
    
    fftw_reverse_plan_ = fftwf_plan_dft_2d(L_, L_, fftw_in_,
                                           fftw_out_,
                                           FFTW_BACKWARD, FFTW_MEASURE);
    
    ArrayXcf yn_scale_vec = ArrayXcf::Zero(L_);
    for(int n = 0; n < N_; n++) {
        yn_scale_vec(n) = std::pow(A_, -n) * std::pow(W_, (n*n)/2.0f); 
    }
    // http://stackoverflow.com/a/20515413/1073963  is noailais() going to save
    // us anything here?
    
    yn_scale_ = yn_scale_vec.matrix() * yn_scale_vec.matrix().transpose(); 
    
    ArrayXcf vn = ArrayXcf::Zero(L_);
    for(int n = 0; n < M_; ++n) {
        vn(n) = std::pow(W, (-n*n)/2.0f);
    }
    
    for(int n = L_-N_+1; n < L_; ++n) {
        vn(n) = std::pow(W_, -((L_-n)*(L_-n))/2.0f); 
    }

    Vr_ = fft(vn.matrix() * vn.matrix().transpose());
            
    ArrayXcf g_scale_vec = ArrayXcf::Zero(M_);
    for(int k = 0; k < M_; ++k) {
        g_scale_vec(k) = std::pow(W_, (k*k/2.0f)); 
    }

    g_scale_ = g_scale_vec.matrix() * g_scale_vec.matrix().transpose(); 

}


ChirpZ2d::~ChirpZ2d(void )
{
    fftwf_free(fftw_in_);
    fftwf_free(fftw_out_); 
    
}

MatrixXcf ChirpZ2d::compute(const MatrixXcf & x) {

    MatrixXcf yn = MatrixXcf::Zero(L_, L_);
    
    for(int i = 0; i < x.rows(); ++i) {
        for(int j = 0; j < x.cols(); ++j ) {
            // FIXME IS THIS THE RIGHT ORDER? memory matters
            yn(i, j) = x(i, j) * yn_scale_(i, j);
        }
    }

    
    MatrixXcf Yr = fft(yn);

    MatrixXcf Gr = (Yr.array() * Vr_.array()).matrix();

    MatrixXcf gk = ifft(Gr);
    
    MatrixXcf Xk = MatrixXcf::Zero(M_, M_);
    for (int i = 0; i < M_; ++i) {
        for(int j = 0; j < M_; ++j) { 
            Xk(i, j) = g_scale_(i, j) * gk(i, j);
        }
    }

    MatrixXcf out =  Xk / (2*N_*2*N_);

    return out; 

}

MatrixXcf ChirpZ2d::fft(const MatrixXcf & in) {
    assert(in.rows() == L_); 
    assert(in.cols() == L_); 
    Map<MatrixXcf> in_f((fc32_t *)fftw_in_, L_, L_);

    in_f = in;

    fftwf_execute(fftw_forward_plan_);
    
    Map<MatrixXcf> out_field((fc32_t*) fftw_out_, L_, L_) ;
    MatrixXcf out = out_field;
    return out; 
}

MatrixXcf ChirpZ2d::ifft(const MatrixXcf & in) {
    assert(in.rows() == L_); 
    assert(in.cols() == L_); 

    Map<MatrixXcf> in_f((fc32_t *)fftw_in_, L_, L_);
    in_f = in;
    fftwf_execute(fftw_reverse_plan_);
    
    Map<MatrixXcf> out_field((fc32_t*) fftw_out_, L_, L_) ;
    MatrixXcf out = out_field;
    return out;
    

}

}

