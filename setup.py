#!/usr/bin/python
"""\

GAT
***

Genomic annotation tool
  
"""\

import os
import sys
import glob
import numpy

from distribute_setup import use_setuptools
use_setuptools()

from setuptools import Extension, setup, find_packages

#######################################################
#######################################################
try:
    from Cython.Distutils import build_ext
except ImportError:
    # no Cython available - use existing C code
    cmdclass = { }
    GatSegmentList_sources = ['GatSegmentList/GatSegmentList.c', 
                              'utils/gat_utils.c']
    GatEngine_sources = ['GatEngine/GatEngine.c', 
                         'utils/gat_utils.c' ]
else:
    cmdclass = { 'build_ext' : build_ext }
    GatSegmentList_sources = ['GatSegmentList/GatSegmentList.pyx', 
                              'utils/gat_utils.c']
    GatEngine_sources = ['GatEngine/GatEngine.pyx', 
                         'utils/gat_utils.c' ]

name = "gat"
version = "0.2"

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

# link against rt for shared memory access
GatSegmentList = Extension(
   "GatSegmentList",                   # name of extension
   GatSegmentList_sources,
   libraries=[ "z", 'rt' ],
   library_dirs = [],
   include_dirs=['./utils', numpy.get_include() ],
   language="c",
   )

GatEngine = Extension(
    "GatEngine",                   # name of extension
    GatEngine_sources,
    libraries=[ "z", 'rt' ],
    library_dirs = [],
    include_dirs=["./utils", 
                  "/usr/lib64/python2.6/site-packages/numpy/core/include", 
                  numpy.get_include() ],
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
    'ext_modules': [ GatSegmentList, GatEngine ],
    'cmdclass' : {'build_ext': build_ext },
    'scripts' : ['scripts/gat-run.py', 
                 'scripts/gat-great.py',
                 'scripts/gat-compare.py',
                 'scripts/gat-plot.py' ],
    'zip_safe' :False,
   }

if __name__=='__main__':
   dist = setup(**metadata)
