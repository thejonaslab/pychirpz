

all: forwardtest foo.o

forwardtest : forwardtest.cc
	g++-5 -std=c++14 -g -O3 -ffast-math -march=native forwardtest.cc  -lm -lfftw3f -lboost_system -lboost_timer -I/usr/local/include/eigen3 -I./eigen -o forwardtest

foo.o : foo.cc
	g++-5 -std=c++14 -g -O3 -ffast-math -march=native foo.cc  -lm -lfftw3f -lboost_system -lboost_timer -I/usr/local/include/eigen3 -I./eigen 
