# distutils: language = c++
# distutils: sources = chirpz.cc

cimport numpy as np

from eigency.core cimport *

cdef extern from "chirpz.h" namespace "chirpz":

     cdef void hello()


def hello_world():
    hello()

