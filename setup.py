#!/usr/bin/python
"""\

GAT
***

Genomic annotation tool

"""\

import sys
import glob
import platform

########################################################################
#######################################################################
# Check for dependencies
##
# Is there a way to do this more elegantly?
# 1. Run "pip install numpy"
# 2. Wrap inside functions (works for numpy/pysam, but not cython)
try:
    import numpy
except ImportError:
    raise ImportError(
        "GAT requires numpy to be installed before running setup.py "
        "(pip install numpy)")

try:
    import Cython
except ImportError:
    raise ImportError(
        "GAT code requires cython to be installed before running setup.py "
        "(pip install cython)")

########################################################################
########################################################################
# Import setuptools
# Use existing setuptools, otherwise try ez_setup.
try:
    import setuptools
except ImportError:
    # try to get via ez_setup
    # ez_setup did not work on all machines tested as
    # it uses curl with https protocol, which is not
    # enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    raise ImportError(
        "GAT requires setuptools 1.1 higher")

###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major == 2 and minor1 < 5 or major < 2:
    raise SystemExit("""GAT requires Python 2.5 or later.""")

try:
    from Cython.Distutils import build_ext
except ImportError:
    # no Cython available - use existing C code
    cmdclass = {}
    CoordinateList_sources = ['gat/CoordinateList.c',
                              'utils/gat_utils.c']
    SegmentList_sources = ['gat/SegmentList.c',
                           'utils/gat_utils.c']
    PositionList_sources = ['gat/PositionList.c',
                            'utils/gat_utils.c']
    Engine_sources = ['gat/Engine.c',
                      'utils/gat_utils.c']
else:
    cmdclass = {'build_ext': build_ext}
    CoordinateList_sources = ['gat/CoordinateList.pyx',
                              'utils/gat_utils.c']
    SegmentList_sources = ['gat/SegmentList.pyx',
                           'utils/gat_utils.c']
    PositionList_sources = ['gat/PositionList.pyx',
                            'utils/gat_utils.c']
    Engine_sources = ['gat/Engine.pyx',
                      'utils/gat_utils.c']


install_requires = []

for requirement in (l.strip() for l in open('requirements.txt')
                    if not l.startswith("#")):
    install_requires.append(requirement)


name = "gat"
version = "1.3.6"

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

libraries = ['z']

# link against rt for shared memory access
# macosx has no librt
if platform.system() in ("Linux", "Windows"):
    libraries.append('rt')

CoordinateList = Extension(
    "gat.CoordinateList",
    CoordinateList_sources,
    libraries=libraries,
    library_dirs=[],
    include_dirs=['./utils', numpy.get_include()],
    language="c",
)

SegmentList = Extension(
    "gat.SegmentList",
    SegmentList_sources,
    libraries=libraries,
    library_dirs=[],
    include_dirs=['./utils', numpy.get_include()],
    language="c",
)

PositionList = Extension(
    "gat.PositionList",
    PositionList_sources,
    libraries=libraries,
    library_dirs=[],
    include_dirs=['./utils', numpy.get_include()],
    language="c",
)

Engine = Extension(
    "gat.Engine",
    Engine_sources,
    libraries=libraries,
    library_dirs=[],
    include_dirs=["./utils",
                  "../GatSegmentList",
                  numpy.get_include()],
    language="c",
)

setup(
    name=name,
    version=version,
    description="GenomicAssocationTester",
    long_description=__doc__,
    author="Andreas Heger",
    author_email="andreas.heger@gmail.com",
    license="MIT",
    platforms="ALL",
    url="https://github.com/AndreasHeger/gat",
    py_modules=[x[:-3] for x in glob.glob('gat/*.py')],
    install_requires=install_requires,
    ext_modules=[CoordinateList, SegmentList, PositionList, Engine],
    cmdclass={'build_ext': build_ext},
    scripts=['scripts/gat-run.py',
             'scripts/gat-great.py',
             'scripts/gat-compare.py',
             'scripts/gat-plot.py'],
    zip_safe=False,
)
