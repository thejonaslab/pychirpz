#include <fftw3.h>
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>
#include <cmath>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <boost/program_options.hpp>

/*
  This is a simple stand-alone C++ benchmark of the C++ 
  library to show how you can use it without python
  and to provide a better framework for microbenchmarking. 
 */

int chirpz_benchmark(int argc, const char *argv[]) {
    options_description desc{"Options"};
    desc.add_options()
        ("help,h", "Help screen")
        ("M", value<int>()->default_value(256), "input point number")
        ("N", value<int>()->default_value(256), "output point number")
        ("ITERS", value<int>()->default_value(1000), "Number of transforms to do ")
        ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    int M = vm["M"].as<int>(); 
    int N = vm["N"].as<int>(); 

    int ITERS = vm["ITERS"].as<int>();
    
    float theta_start = -PI/2.0; 
    float theta_delta = PI / 256; 
    
    fc32_t A = std::exp(J * theta_start);
    fc32_t W = std::exp(-J * theta_delta);

    boost::timer::cpu_timer setup_timer;
    chirpz::ChirpZ2d<chirpz::c32_t> cz(N, M, A, W) ;
    auto setup_ns = boost::chrono::nanoseconds(setup_timer.elapsed().wall); 
    auto setup_s = boost::chrono::duration_cast<boost::chrono::seconds>(setup_ns);
    std::cout << setup_s << " for setup" << std::endl;

    {
        MatrixXcf x = MatrixXcf::Random(M, M);
        boost::timer::cpu_timer timer;
        for(int i = 0; i < ITERS; ++i) { 
            auto y = cz.compute(x); 
        }
        auto ns = boost::chrono::nanoseconds(timer.elapsed().wall);
        auto ms = boost::chrono::duration_cast<boost::chrono::milliseconds>(ns);
 
        std::cout << ms / double(ITERS) << " per eval" << std::endl; 
    }

}


int main(int argc, const char *argv[]) {
    //spherical_bench(); 
    //fft_benchmark(argc, argv);

    chirpz_benchmark(argc, argv); 

}
