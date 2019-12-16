#! /bin/bash

set -e

# Retrieve the latest miniconda distribution for linux
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;

# Install miniconda
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

# Configure conda
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda config --add channels conda-forge
conda create -q -n travis_env python=$TRAVIS_PYTHON_VERSION
source activate travis_env

# Install
conda install --file $TRAVIS_BUILD_DIR/requirements.txt
python $TRAVIS_BUILD_DIR/setup.py install
