#!/usr/bin/python
"""\

FastGTF
*******

Module for fast parsing of `gtf <http://mblab.wustl.edu/GTF2.html>`_ formatted files.

  
"""\

import os
import sys
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

name = "gat"
version = "0.1"

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
      library_dirs=[],
      libraries=[],
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
      "gat/__init__", ],
    'ext_modules': [cgat],
    'cmdclass' : {'build_ext': build_ext},
   }

if __name__=='__main__':
   dist = setup(**metadata)
