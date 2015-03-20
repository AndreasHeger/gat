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

   pip install gat

To install on OS X, we suggest to begin by installing 
`homebrew <http://brew.sh/>`_ by following these
`instructions <http://hackercodex.com/guide/mac-osx-mountain-lion-10.8-configuration/>`_

Follow then by::

   brew install python --with-brewed-openssl
   pip install numpy
   pip install cython
   pip install gat   

Installing via source
=====================

The latest changes can be obtained by cloning the repository on github_::

   git clone https://github.com/AndreasHeger/gat.git

To install, type::

   python setup.py install

in the package directory.


Release History
===============

1.2.2 Minor features
   * Added --random-seed as option.
   * moved documentation to read-the-docs.

1.2.1 Bugfix release:
   * added missing files :file:`requires.txt` to tarball

1.2 Bugfix release:
   * Command line options renamed for CGAT compatibility
   * minor bugfixes

1.1 Bugfix release: easier Galaxy integration
   * Changed to distutils (from distribute)
   * Changed /bin/env to /usr/bin/env

1.0 Release coinciding with publication

0.2
  * First release

0.1 
   * Alpha release
