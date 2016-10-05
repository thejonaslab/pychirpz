# distutils: language = c++
# distutils: sources = chirpz.cc

cimport numpy as np
import numpy as np


from eigency.core cimport *

from libcpp.complex cimport complex as cc

cdef extern from "chirpz.h" namespace "chirpz":
     ctypedef fc32_t

     cdef void hello(cc[float])

     cdef cppclass ChirpZ:
         ChirpZ(int, int, cc[float], cc[float])
         ArrayXcf compute(Map[ArrayXcf] &)

     cdef cppclass ChirpZ2d:
         ChirpZ2d(int, int, cc[float], cc[float])
         MatrixXcf compute(Map[MatrixXcf] &)

def hello_world(c):
    hello(c)

cdef class PyChirpZ:
   cdef ChirpZ *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new ChirpZ(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[ArrayXcf](x)))

cdef class PyChirpZ2d:
   cdef ChirpZ2d *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new ChirpZ2d(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcf](x)))

def zoom_fft(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = PyChirpZ(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex64)).flatten()



def zoom_fft2(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = PyChirpZ2d(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex64))

