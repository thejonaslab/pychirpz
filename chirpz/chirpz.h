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
    typedef float number_t ; 
    typedef std::complex<number_t> c_t;
    typedef fftwf_complex fft_complex;
    typedef fftwf_plan fft_plan;
    typedef MatrixXcf Matrix;
    typedef ArrayXcf Array;

    inline static fft_complex * fft_alloc(std::size_t s) {
        return fftwf_alloc_complex(s); 
    }
    
    inline static void  fft_free (void * p) {
        fftwf_free(p);
    }

    inline static void fft_execute(fft_plan p) {
        fftwf_execute(p);
    }
    
    inline static fft_plan fft_plan_dft_1d(int a, fft_complex* b, fft_complex* c,
                             int d, unsigned int e) {
        return fftwf_plan_dft_1d(a, b, c, d, e);
    }

    inline static fft_plan fft_plan_dft_2d(int a, int b, fft_complex* c, fft_complex* d,
                             int e, unsigned int f){
        
        return fftwf_plan_dft_2d(a, b, c, d, e, f);
    }

};

struct c64_t
{ 
    typedef double number_t ; 
    typedef std::complex<number_t> c_t;
    typedef fftw_complex fft_complex;
    typedef fftw_plan fft_plan;
    typedef MatrixXcd Matrix;
    typedef ArrayXcd Array;
    
    inline static fft_complex * fft_alloc(std::size_t s) {
        return fftw_alloc_complex(s); 
    }
    
    inline static void  fft_free (void * p) {
        fftw_free(p);
    }

    inline static void fft_execute(fft_plan p) {
        fftw_execute(p);
    }
    
    inline static fft_plan fft_plan_dft_1d(int a, fft_complex* b,
                                           fft_complex* c,
                                           int d, unsigned int e) {
        return fftw_plan_dft_1d(a, b, c, d, e);
    }

    inline static fft_plan fft_plan_dft_2d(int a, int b, fft_complex* c,
                                           fft_complex* d,
                                           int e, unsigned int f){
        
        return fftw_plan_dft_2d(a, b, c, d, e, f);
    }

};


template <class T>
class ChirpZ  {
    
public:
    typedef typename T::number_t number_t;
    typedef typename T::c_t c_t;
    typedef typename T::Array Array;
    
    ChirpZ(int N, int M, c_t A, c_t W, bool fftw_patient=true) :
        N_(N), M_(M), A_(A), W_(W),
        L_(std::pow(2, nextpo2(N + M -1)))
    {

        // allocate FFT 
        fft_in_ = (typename T::fft_complex*)T::fft_alloc(L_);
        fft_out_ = (typename T::fft_complex*)T::fft_alloc(L_);
        int fftw_measure = FFTW_MEASURE;
        if(fftw_patient) {
            fftw_measure = FFTW_PATIENT;
        }
        
        fft_forward_plan_ = T::fft_plan_dft_1d(L_, fft_in_,
                                               fft_out_,
                                               FFTW_FORWARD, fftw_measure);
        
        fft_reverse_plan_ = T::fft_plan_dft_1d(L_, fft_in_,
                                               fft_out_,
                                               FFTW_BACKWARD, fftw_measure);
        
        yn_scale_ = Array::Zero(L_);
        for(int n = 0; n < N_; n++) {
            c_t a = c_t(std::pow(A_, -n)); 
            c_t b = std::pow(W_, static_cast<number_t>((n*n)/2.0));

            yn_scale_(n) = a * b;
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
    typedef typename T::number_t number_t;
    typedef typename T::Array Array;
    typedef typename T::Matrix Matrix;

public:
    ChirpZ2d(int N, int M, c_t A, c_t W, bool fftw_patient=true) :
        N_(N), M_(M), A_(A), W_(W),
        L_(std::pow(2, nextpo2(N + M -1)))
    {
        
        // allocate FFT 
        fft_in_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);
        fft_out_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);
        ifft_in_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);
        ifft_out_ = (typename T::fft_complex*)T::fft_alloc(L_ * L_);

        int fftw_measure = FFTW_MEASURE;
        if(fftw_patient) {
            fftw_measure = FFTW_PATIENT;
        }
        fft_forward_plan_ = T::fft_plan_dft_2d(L_, L_, fft_in_,
                                               fft_out_,
                                               FFTW_FORWARD, fftw_measure); 
        
        fft_reverse_plan_ = T::fft_plan_dft_2d(L_, L_, ifft_in_,
                                               ifft_out_,
                                               FFTW_BACKWARD, fftw_measure);
        Array yn_scale_vec = Array::Zero(L_);
        for(int n = 0; n < N_; n++) {
            c_t a = c_t(std::pow(A_, -n)); 
            c_t b = std::pow(W_, static_cast<number_t>((n*n)/2.0));

            yn_scale_vec(n) = a * b; 
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
        
        g_scale_ = g_scale_vec.matrix() * g_scale_vec.matrix().transpose() / (2*N_*2*N_); 
        
    }
    
    
    ~ChirpZ2d(void )
    {
        T::fft_free(fft_in_);
        T::fft_free(fft_out_); 
        T::fft_free(ifft_in_);
        T::fft_free(ifft_out_); 
        
    }
    
    Matrix compute(const Matrix & x) {


        Map<Matrix> in_f((c_t *)fft_in_, L_, L_);
        in_f = Matrix::Zero(L_, L_);
        in_f.block(0, 0, x.rows(), x.cols()) =
            (x.array() * yn_scale_.block(0, 0, x.rows(), x.cols()).array()).matrix();
        
        
        T::fft_execute(fft_forward_plan_);
        
        Map<Matrix> out_field((c_t*) fft_out_, L_, L_) ;

        Map<Matrix> ifft_in_f((c_t *)ifft_in_, L_, L_);
        ifft_in_f = (out_field.array() * Vr_.array()).matrix();
        
        T::fft_execute(fft_reverse_plan_);
        
        Map<Matrix> gk((c_t*) ifft_out_, L_, L_) ;
        
        auto Xk = (g_scale_.array() * gk.block(0, 0, M_, M_).array()).matrix();
        
        return Xk; 
        
    }

private:
    
    Matrix fft(const Ref<const Matrix> & in) {
        assert(in.rows() == L_); 
        assert(in.cols() == L_); 
        Map<Matrix> in_f((c_t *)fft_in_, L_, L_);
        
        in_f = in;
        
        T::fft_execute(fft_forward_plan_);
        
        Map<Matrix> out_field((c_t*) fft_out_, L_, L_) ;
        //Matrix out = out_field;
        return out_field; 
    }
    
    Matrix ifft(const Ref<const Matrix> & in) {
        assert(in.rows() == L_); 
        assert(in.cols() == L_); 
        
        Map<Matrix> in_f((c_t *)ifft_in_, L_, L_);
        in_f = in;
        T::fft_execute(fft_reverse_plan_);
        
        Map<Matrix> out_field((c_t*) ifft_out_, L_, L_) ;
        //Matrix out = out_field;
        return out_field;
        
        
    }
    

    int N_; 
    int M_;
    c_t A_;
    c_t W_;

    int L_;
    Matrix yn_scale_; 

    typename T::fft_complex* fft_in_;
    typename T::fft_complex* fft_out_;
    
    typename T::fft_complex* ifft_in_;
    typename T::fft_complex* ifft_out_;
    
    typename T::fft_plan fft_forward_plan_;
    typename T::fft_plan fft_reverse_plan_;

    Matrix Vr_;
    Matrix g_scale_;

        

}; 

    
}
