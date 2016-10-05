#!/usr/bin/env python
import numpy as np
import os
import eigency
from distutils.core import setup
#import pkgconfig

from distutils.extension import Extension
from Cython.Build import cythonize

#os.environ["CC"] = "g++-5" 
#os.environ["CXX"] = "g++-5"

eigen_include = "" # pkgconfig.cflags('eigen3')
    
extensions = [Extension(name = "chirpz.cychirpz",
                        sources = ["chirpz/cychirpz.pyx", "chirpz/chirpz.cc"],
                        extra_compile_args = ['-O3', '-DEIGEN_NO_AUTOMATIC_RESIZING', '-march=native', '-std=c++14', eigen_include],
                        include_dirs = [np.get_include(), "./",] \
                        + eigency.get_includes(include_eigen=True), 
                        extra_link_args = ['-lm', '-lfftw3f', '-lboost_system', '-lboost_timer'], 
                        language='c++')
]


setup(name='chirpz',
      version='1.0',
      description='Python and C++ implementation of Chirp-Z transform',
      author='Eric Jonas',
      author_email='jonas@eecs.berkeley.edu',
      url='https://www.github.com/ericmjonas/pychirpz/',
      packages=['chirpz'],
      ext_modules = cythonize(extensions),
     )
