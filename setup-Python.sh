#!/usr/bin/env bash

# Create virtual environment
pip install virtualenv
virtualenv cgat-venv
source cgat-venv/bin/activate

# Install some Python prerequisites
pip install Cython
pip install numpy
pip install matplotlib
pip install scipy
pip install -r https://raw.github.com/AndreasHeger/sphinx-report/master/requires.txt
pip install --upgrade setuptools


