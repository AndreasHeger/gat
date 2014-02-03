#!/usr/bin/env bash

# Build Python
cd
mkdir CGAT
wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
tar xzvf Python-2.7.5.tgz
rm Python-2.7.5.tgz
cd Python-2.7.5
./configure --prefix=$HOME/CGAT/Python-2.7.5
make
make install
cd
rm -rf Python-2.7.5

# Create virtual environment
cd
cd CGAT
wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
tar xvfz virtualenv-1.10.1.tar.gz
rm virtualenv-1.10.1.tar.gz
cd virtualenv-1.10.1
$HOME/CGAT/Python-2.7.5/bin/python virtualenv.py cgat-venv
source cgat-venv/bin/activate

# Install some Python prerequisites
pip install Cython
pip install numpy
pip install matplotlib
pip install scipy
pip install -r https://raw.github.com/AndreasHeger/sphinx-report/master/requires.txt
pip install --upgrade setuptools


