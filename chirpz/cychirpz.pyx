# distutils: language = c++

cimport numpy as np
import numpy as np


from eigency.core cimport *

from libcpp.complex cimport complex as cc

cdef extern from "chirpz.h" namespace "chirpz":



     cdef void hello(cc[float])

     cdef struct c32_t:
         pass

     cdef struct c64_t:
        pass

     cdef cppclass ChirpZ32 "chirpz::ChirpZ<chirpz::c32_t>":
         ChirpZ32(int, int, cc[float], cc[float])
         ArrayXcf compute(Map[ArrayXcf] &)


     cdef cppclass ChirpZ64 "chirpz::ChirpZ<chirpz::c64_t>":
         ChirpZ64(int, int, cc[double], cc[double])
         ArrayXcd compute(Map[ArrayXcd] &)



     cdef cppclass ChirpZ2d32 "chirpz::ChirpZ2d<chirpz::c32_t>":
         ChirpZ2d32(int, int, cc[float], cc[float])
         MatrixXcf compute(Map[MatrixXcf] &)


     cdef cppclass ChirpZ2d64 "chirpz::ChirpZ2d<chirpz::c64_t>":
         ChirpZ2d64(int, int, cc[double], cc[double])
         MatrixXcd compute(Map[MatrixXcd] &)


def hello_world(c):
    hello(c)

cdef class PyChirpZ32:
   cdef ChirpZ32 *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new ChirpZ32(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[ArrayXcf](x)))


cdef class PyChirpZ64:
   cdef ChirpZ64 *thisptr

   def __cinit__(self, int N, int M, cc[double] A, cc[double] W):
       self.thisptr = new ChirpZ64(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[ArrayXcd](x)))


cdef class PyChirpZ2d32:
   cdef ChirpZ2d32 *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new ChirpZ2d32(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcf](x)))

cdef class PyChirpZ2d64:
   cdef ChirpZ2d64 *thisptr

   def __cinit__(self, int N, int M, cc[double] A, cc[double] W):
       self.thisptr = new ChirpZ2d64(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcd](x)))



def zoom_fft(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = PyChirpZ64(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex128)).flatten()



def zoom_fft2(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = PyChirpZ2d64(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex128))

