# distutils: language = c++

cimport numpy as np
import numpy as np


from eigency.core cimport *

from libcpp.complex cimport complex as cc

cdef extern from "chirpz.h" namespace "chirpz":



     cdef void hello(cc[float])

     cdef cppclass ChirpZ:
         ChirpZ(int, int, cc[float], cc[float])
         ArrayXcf compute(Map[ArrayXcf] &)

     cdef struct c32_t:
         pass

     cdef struct c64_t:
        pass

     cdef cppclass TChirpZ32 "chirpz::TChirpZ<chirpz::c32_t>":
         TChirpZ32(int, int, cc[float], cc[float])
         ArrayXcf compute(Map[ArrayXcf] &)


     cdef cppclass TChirpZ64 "chirpz::TChirpZ<chirpz::c64_t>":
         TChirpZ64(int, int, cc[double], cc[double])
         ArrayXcd compute(Map[ArrayXcd] &)


     cdef cppclass ChirpZ2d:
         ChirpZ2d(int, int, cc[float], cc[float])
         MatrixXcf compute(Map[MatrixXcf] &)

     cdef cppclass TChirpZ2d32 "chirpz::TChirpZ2d<chirpz::c32_t>":
         TChirpZ2d32(int, int, cc[float], cc[float])
         MatrixXcf compute(Map[MatrixXcf] &)


     cdef cppclass TChirpZ2d64 "chirpz::TChirpZ2d<chirpz::c64_t>":
         TChirpZ2d64(int, int, cc[double], cc[double])
         MatrixXcd compute(Map[MatrixXcd] &)


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

cdef class TPyChirpZ32:
   cdef TChirpZ32 *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new TChirpZ32(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[ArrayXcf](x)))


cdef class TPyChirpZ64:
   cdef TChirpZ64 *thisptr

   def __cinit__(self, int N, int M, cc[double] A, cc[double] W):
       self.thisptr = new TChirpZ64(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[ArrayXcd](x)))


cdef class PyChirpZ2d:
   cdef ChirpZ2d *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new ChirpZ2d(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcf](x)))

cdef class TPyChirpZ2d32:
   cdef TChirpZ2d32 *thisptr

   def __cinit__(self, int N, int M, cc[float] A, cc[float] W):
       self.thisptr = new TChirpZ2d32(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcf](x)))

cdef class TPyChirpZ2d64:
   cdef TChirpZ2d64 *thisptr

   def __cinit__(self, int N, int M, cc[double] A, cc[double] W):
       self.thisptr = new TChirpZ2d64(N, M, A, W)
       
   def __dealloc__(self):
       del self.thisptr

   def compute(self, x):
       return ndarray(self.thisptr.compute(Map[MatrixXcd](x)))



def zoom_fft(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = TPyChirpZ64(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex128)).flatten()



def zoom_fft2(x, theta_start, step_size, M):
    
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    pCZ = TPyChirpZ2d64(len(x), M, A, W)
    
    return pCZ.compute(x.astype(np.complex128))

