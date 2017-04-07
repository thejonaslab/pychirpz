#!/bin/bash

set -ev

brew update
brew install fftw boost

wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;

exit 0;
