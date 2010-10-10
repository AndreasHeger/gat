#!/usr/bin/python
"""\

GAT
***

Genomic annotation tool
  
"""\

import os
import sys
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

name = "gat"
version = "0.1.2"

classifiers = """
Development Status :: 5 - Production/Stable
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

cgat = Extension(
    "cgat",                   # name of extension
    [ "gat/cgat.pyx", "gat/gat_utils.c" ],
      libraries=[ "z" ],
      include_dirs=["/cpp-software/lib/python2.6/site-packages/numpy/core/include",],
      language="c",
    )


metadata = {
    'name': name,
    'version': version,
    'description': "GenomicEnrichmentTool", 
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "GPL",
    'platforms': "ALL",
    'url': "",
    'py_modules': [
      "gat/__init__", "gat/Stats", "gat/IOTools", "gat/Experiment", "gat/Bed" ],
    'ext_modules': [cgat],
    'cmdclass' : {'build_ext': build_ext},
    'scripts' : ['scripts/gatrun.py' ],
   }

if __name__=='__main__':
   dist = setup(**metadata)
