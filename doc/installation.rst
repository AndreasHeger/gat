============
Installation
============

Requirements
============

GAT has been written in `python <http://www.python.org>`_ and has been
tested with the following python versions:

   * python 2.7.3
   * python 2.6.8

GAT requires the following modules to be installed at installation:

   * `numpy <http://www.numpy.org/>`_ 1.4 or greater
   * `cython <http://www.cython.org/>`_ 0.14 or greater

The plotting and unit test modules also require scipy and matplotlib.

Installing from PyPi
====================

GAT is available at the `python package index
<https://pypi.python.org/pypi>`_ and can be installed
using `pip <http://www.pip-installer.org/en/latest/>`_ or 
`setuptools <https://pypi.python.org/pypi/setuptools>`_.

To install via pip, type::

   pip install genomic-association-tester

Installing via source
=====================

The latest released version of Gat can be downloaded from the following URL:

   http://code.google.com/p/genomic-association-tester/downloads/list

Next, uncompress the tar archive::
   
   tar -xvzf genomic-association-tester.tar.gz

The latest changes can be obtained by cloning the repository
using merucial::

   hg clone https://andreas.heger@code.google.com/p/genomic-association-tester/ gat

To install, type::

   python setup.py

in the package directory.


Release History
===============

0.2
  * First release

0.1 
   * Alpha release
