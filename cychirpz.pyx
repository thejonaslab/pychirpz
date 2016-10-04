# distutils: language = c++
# distutils: sources = chirpz.cc

cimport numpy as np

from eigency.core cimport *

from libcpp.complex cimport complex as cc

cdef extern from "chirpz.h" namespace "chirpz":
     ctypedef fc32_t

     cdef void hello(cc[float])

     cdef cppclass ChirpZ:
         ChirpZ(int, int, cc[float], cc[float])
         ArrayXcf compute(Map[ArrayXcf] &)

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
