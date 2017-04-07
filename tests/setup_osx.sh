#!/bin/bash

set -ev

brew update
brew install fftw || true
# if it's already installed it will fail, so we always return success
brew install boost || true 

wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;

exit 0;
