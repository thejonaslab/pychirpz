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
inline int nextpo2(int x) {
    return int(std::ceil(std::log2(x)));
    
}

// use a traits class
struct c32_t
{ 
    typedef std::complex<float> c_t;
    typedef fftwf_complex fft_complex;
    typedef fftwf_plan fft_plan;
    typedef MatrixXcf Matrix;
    typedef ArrayXcf Array;
    static constexpr fftwf_complex* (*fft_alloc) (std::size_t)  = fftwf_alloc_complex;
    static constexpr void (*fft_free) (void *) = fftwf_free;
    static constexpr void (*fft_execute) (fft_plan) = fftwf_execute;
    static constexpr fft_plan (*fft_plan_dft_1d)(int, fft_complex*, fft_complex*,
                                                 int, unsigned int)  = fftwf_plan_dft_1d; 

    static constexpr fft_plan (*fft_plan_dft_2d)(int, int, fft_complex*, fft_complex*,
                                                 int, unsigned int)  = fftwf_plan_dft_2d; 

};

struct c64_t
{ 
    typedef std::complex<double> c_t;
    typedef fftw_complex fft_complex;
    typedef fftw_plan fft_plan;
    typedef MatrixXcd Matrix;
    typedef ArrayXcd Array;
    static constexpr fftw_complex* (*fft_alloc) (std::size_t)  = fftw_alloc_complex;
    static constexpr void (*fft_free) (void *) = fftw_free;
    static constexpr void (*fft_execute) (fft_plan) = fftw_execute;
    static constexpr fft_plan (*fft_plan_dft_1d)(int, fft_complex*, fft_complex*,
                                                 int, unsigned int)  = fftw_plan_dft_1d; 

    static constexpr fft_plan (*fft_plan_dft_2d)(int, int, fft_complex*, fft_complex*,
                                                 int, unsigned int)  = fftw_plan_dft_2d; 
};


template <class T>
class ChirpZ  {
    
public:
    typedef typename T::c_t c_t;
    typedef typename T::Array Array;
    
    ChirpZ(int N, int M, c_t A, c_t W) :
        N_(N), M_(M), A_(A), W_(W),
        L_(std::pow(2, nextpo2(N + M -1)))
    {

        // allocate FFT 
        fft_in_ = (typename T::fft_complex*)T::fft_alloc(L_);
        fft_out_ = (typename T::fft_complex*)T::fft_alloc(L_);
        
        fft_forward_plan_ = T::fft_plan_dft_1d(L_, fft_in_,
                                               fft_out_,
                                               FFTW_FORWARD, FFTW_MEASURE);
        
        fft_reverse_plan_ = T::fft_plan_dft_1d(L_, fft_in_,
                                               fft_out_,
                                               FFTW_BACKWARD, FFTW_MEASURE);
        
        yn_scale_ = Array::Zero(L_);
        for(int n = 0; n < N_; n++) {
            yn_scale_(n) = std::pow(A_, -n) * std::pow(W_, (n*n)/2.0f); 
        }
        
        Array vn = Array::Zero(L_);
        for(int n = 0; n < M_; ++n) {
            vn(n) = std::pow(W, (-n*n)/2.0f);
        }
        
        for(int n = L_-N_+1; n < L_; ++n) {
            vn(n) = std::pow(W_, -((L_-n)*(L_-n))/2.0f); 
        }
        
        Vr_ = fft(vn);
        
        g_scale_ = Array::Zero(M_);
        for(int k = 0; k < M_; ++k) {
            g_scale_(k) = std::pow(W_, (k*k/2.0f)); 
        }
        
        
    
    }
    
    
    ~ChirpZ(void )
    {
        T::fft_free(fft_in_);
        T::fft_free(fft_out_); 
        
        
    }
    
    Array compute(const Array & x) {
        
    
        Array yn = Array::Zero(L_);
        
        for(int i = 0; i < x.size(); ++i) { 
            yn(i) = x(i) * yn_scale_(i); 
        }
        
        
        Array Yr = fft(yn);
        
        Array Gr = Yr * Vr_;
        
        
        Array gk = ifft(Gr);
        
        Array Xk = Array::Zero(g_scale_.size());
        for (int i = 0; i < g_scale_.size(); ++i) {
            Xk(i) = g_scale_(i) * gk(i);
        }
        Array out =  Xk / (2*N_);
        
        return out; 
        
    }
    
    
    
private:
    int N_; 
    int M_;
    c_t A_;
    c_t W_;

    int L_;
    Array yn_scale_; 

    typename T::fft_complex* fft_in_;
    typename T::fft_complex* fft_out_;
    
    typename T::fft_plan fft_forward_plan_;
    typename T::fft_plan fft_reverse_plan_;

    Array Vr_;
    Array g_scale_;

    // convenience functions, warning, make copies
    Array fft(const Array & in) {
        Map<Array> in_f((c_t *)fft_in_, L_);
        
        in_f = in;
        
        T::fft_execute(fft_forward_plan_);
        
        Map<Array> out_field((c_t*) fft_out_, L_) ;
        Array out = out_field;
        return out; 
    }
    
    Array ifft(Array in) {
        
        Map<Array> in_f((c_t *)fft_in_, L_);
        in_f = in;
        T::fft_execute(fft_reverse_plan_);
        
        Map<Array> out_field((c_t*) fft_out_, L_) ;
        Array out = out_field;
        return out;

    }        

}; 

template <class T>
class ChirpZ2d  {
    
public:
    typedef typename T::c_t c_t;
    typedef typename T::Array Array;
    typedef typename T::Matrix Matrix;
    
public:
    ChirpZ2d(int N, int M, c_t A, c_t W) :
        N_(N), M_(M), A_(A), W_(W),
        L_(std::pow(2, nextpo2(N + M -1)))
    {
        
        // allocate FFT 
        fft_in_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);
        fft_out_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);
        
        fft_forward_plan_ = T::fft_plan_dft_2d(L_, L_, fft_in_,
                                                fft_out_,
                                                FFTW_FORWARD, FFTW_MEASURE);
        
        fft_reverse_plan_ = T::fft_plan_dft_2d(L_, L_, fft_in_,
                                               fft_out_,
                                               FFTW_BACKWARD, FFTW_MEASURE);
        
        Array yn_scale_vec = Array::Zero(L_);
        for(int n = 0; n < N_; n++) {
            yn_scale_vec(n) = std::pow(A_, -n) * std::pow(W_, (n*n)/2.0f); 
        }
        // http://stackoverflow.com/a/20515413/1073963  is noailais() going to save
        // us anything here?
        
        yn_scale_ = yn_scale_vec.matrix() * yn_scale_vec.matrix().transpose(); 
        
        Array vn = Array::Zero(L_);
        for(int n = 0; n < M_; ++n) {
            vn(n) = std::pow(W, (-n*n)/2.0f);
        }
        
        for(int n = L_-N_+1; n < L_; ++n) {
            vn(n) = std::pow(W_, -((L_-n)*(L_-n))/2.0f); 
        }
        
        Vr_ = fft(vn.matrix() * vn.matrix().transpose());
        
        Array g_scale_vec = Array::Zero(M_);
        for(int k = 0; k < M_; ++k) {
            g_scale_vec(k) = std::pow(W_, (k*k/2.0f)); 
        }
        
        g_scale_ = g_scale_vec.matrix() * g_scale_vec.matrix().transpose(); 
        
    }
    
    
    ~ChirpZ2d(void )
    {
        T::fft_free(fft_in_);
        T::fft_free(fft_out_); 
        
    }
    
    Matrix compute(const Matrix & x) {

        Matrix yn = Matrix::Zero(L_, L_);
        
        for(int i = 0; i < x.rows(); ++i) {
            for(int j = 0; j < x.cols(); ++j ) {
                // FIXME IS THIS THE RIGHT ORDER? memory matters
                yn(i, j) = x(i, j) * yn_scale_(i, j);
            }
        }
        
        
        Matrix Yr = fft(yn);
        
        Matrix Gr = (Yr.array() * Vr_.array()).matrix();
        
        Matrix gk = ifft(Gr);
        
        Matrix Xk = Matrix::Zero(M_, M_);
        for (int i = 0; i < M_; ++i) {
            for(int j = 0; j < M_; ++j) { 
                Xk(i, j) = g_scale_(i, j) * gk(i, j);
            }
        }
        
        Matrix out =  Xk / (2*N_*2*N_);
        
        return out; 
        
    }

private:
    
    Matrix fft(const Matrix & in) {
        assert(in.rows() == L_); 
        assert(in.cols() == L_); 
        Map<Matrix> in_f((c_t *)fft_in_, L_, L_);
        
        in_f = in;
        
        T::fft_execute(fft_forward_plan_);
        
        Map<Matrix> out_field((c_t*) fft_out_, L_, L_) ;
        Matrix out = out_field;
        return out; 
    }
    
    Matrix ifft(const Matrix & in) {
        assert(in.rows() == L_); 
        assert(in.cols() == L_); 
        
        Map<Matrix> in_f((c_t *)fft_in_, L_, L_);
        in_f = in;
        T::fft_execute(fft_reverse_plan_);
        
        Map<Matrix> out_field((c_t*) fft_out_, L_, L_) ;
        Matrix out = out_field;
        return out;
        
        
    }
    

    int N_; 
    int M_;
    c_t A_;
    c_t W_;

    int L_;
    Matrix yn_scale_; 

    typename T::fft_complex* fft_in_;
    typename T::fft_complex* fft_out_;
    
    typename T::fft_plan fft_forward_plan_;
    typename T::fft_plan fft_reverse_plan_;

    Matrix Vr_;
    Matrix g_scale_;

        

}; 

    
}
