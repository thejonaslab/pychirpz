

all: cppbench

CHIRPZ = ../chirpz/
EIGEN=/data/jonas/anaconda/lib/python2.7/site-packages/eigency/eigen_3.2.8/

CXX=g++
CXXFLAGS = -std=c++14 -g -O3 -ffast-math -march=native -I${EIGEN} -I ${CHIRPZ} -I./eigen
LDFLAGS= -lfftw3f -lfftw3 -lboost_system -lboost_timer -lboost_program_options -lm

cppbench.o : cppbench.cc
	${CXX} -c cppbench.cc  -o cppbench.o ${CXXFLAGS}


cppbench: cppbench.o
	${CXX} cppbench.o -o cppbench ${CXXFLAGS} ${LDFLAGS}  

clean:
	rm cppbench.o
