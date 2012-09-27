#!/usr/bin/python
"""\

GAT
***

Genomic annotation tool
  
"""\

import os
import sys
import glob

#from ez_setup import use_setuptools
#use_setuptools()

#from setuptools import Extension, setup

from distutils.core import setup, Extension
from Cython.Distutils import build_ext

## note that for automatic cythoning, 
## both pyrex and cython need to be installed.
## see http://mail.python.org/pipermail/distutils-sig/2007-September/008204.html
#try:
#    from Cython.Distutils import build_ext
#except:
#    from setuptools.command.build_ext import build_ext

name = "gat"
version = "0.1"

classifiers = """
Development Status :: 4 - Beta
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""

segmentlist = Extension(
    "csegmentlist",                   # name of extension
    [ "segmentlist/csegmentlist.pyx", "utils/gat_utils.c" ],
    libraries=[ "z", ],
    library_dirs = [],
    include_dirs=['./utils', "/usr/lib64/python2.6/site-packages/numpy/core/include", '.'],
    language="c",
    )

cgat = Extension(
    "cgat",                   # name of extension
    [ "engine/cgat.pyx", "utils/gat_utils.c" ],
    libraries=[ "z" ],
    library_dirs = [],
    include_dirs=["./utils", "/usr/lib64/python2.6/site-packages/numpy/core/include", 
                  'segmentlist', '.'],
    language="c",
    )

metadata = {
    'name': name,
    'version': version,
    'description': "GenomicAssocationTester", 
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "MIT",
    'platforms': "ALL",
    'url': "http://code.google.com/p/genomic-association-tester/",
    'py_modules': [ x[:-3] for x in glob.glob( 'gat/*.py') ],
    'requires' : ['cython (>=0.12)'],
    'ext_modules': [ segmentlist, cgat ],
    'cmdclass' : {'build_ext': build_ext },
    'scripts' : ['scripts/gat-run.py', 
                 'scripts/gat-great.py',
                 'scripts/gat-plot.py' ],
    'zip_safe' :False,

   }

if __name__=='__main__':
   dist = setup(**metadata)
