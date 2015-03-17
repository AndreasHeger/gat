#!/usr/bin/python
"""\

GAT
***

Genomic annotation tool
  
"""\

import os
import sys
import glob
import platform
import re

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
        "GAT requires numpy to be installed before running setup.py (pip install numpy)")

try:
    import Cython
except ImportError:
    raise ImportError(
        "GAT code requires cython to be installed before running setup.py (pip install cython)")

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
    ## enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion, StrictVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")

from Cython.Distutils import build_ext

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

#####################################################################
#####################################################################
# Code to install dependencies from a repository
#####################################################################
# Modified from http://stackoverflow.com/a/9125399
#####################################################################


def which(program):
    """
    Detect whether or not a program is installed.
    Thanks to http://stackoverflow.com/a/377028/70191
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

REPO_REQUIREMENT = re.compile(
    r'^-e (?P<link>(?P<vcs>git|svn|hg|bzr).+#egg=(?P<package>.+)-(?P<version>\d(?:\.\d)*))$')
HTTPS_REQUIREMENT = re.compile(
    r'^-e (?P<link>.*).+#(?P<package>.+)-(?P<version>\d(?:\.\d)*)$')
install_requires = []
dependency_links = []

for requirement in (l.strip() for l in open('requires.txt') if not l.startswith("#")):
    match = REPO_REQUIREMENT.match(requirement)
    if match:
        assert which(match.group('vcs')) is not None, \
            "VCS '%(vcs)s' must be installed in order to install %(link)s" % match.groupdict(
        )
        install_requires.append("%(package)s==%(version)s" % match.groupdict())
        dependency_links.append(match.group('link'))
        continue

    if requirement.startswith("https"):
        install_requires.append(requirement)
        continue

    match = HTTPS_REQUIREMENT.match(requirement)
    if match:
        install_requires.append("%(package)s>=%(version)s" % match.groupdict())
        dependency_links.append(match.group('link'))
        continue

    install_requires.append(requirement)

if major == 2 and minor1 < 5 or major < 2:
    raise SystemExit("""SphinxReport requires Python 2.5 or later.""")

#######################################################
#######################################################
try:
    from Cython.Distutils import build_ext
except ImportError:
    # no Cython available - use existing C code
    cmdclass = {}
    GatSegmentList_sources = ['GatSegmentList/GatSegmentList.c',
                              'utils/gat_utils.c']
    GatEngine_sources = ['GatEngine/GatEngine.c',
                         'utils/gat_utils.c']
else:
    cmdclass = {'build_ext': build_ext}
    GatSegmentList_sources = ['GatSegmentList/GatSegmentList.pyx',
                              'utils/gat_utils.c']
    GatEngine_sources = ['GatEngine/GatEngine.pyx',
                         'utils/gat_utils.c']

name = "gat"
version = "1.2.1"

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


#######################################################################

libraries = ['z']

# macosx has no librt
if platform.system() in ("Linux", "Windows"):
    libraries.append('rt')

# link against rt for shared memory access
GatSegmentList = Extension(
    "GatSegmentList",                   # name of extension
    GatSegmentList_sources,
    libraries=libraries,
    library_dirs=[],
    include_dirs=['./utils', numpy.get_include()],
    language="c",
)

GatEngine = Extension(
    "GatEngine",                   # name of extension
    GatEngine_sources,
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
    ext_modules=[GatSegmentList, GatEngine],
    cmdclass={'build_ext': build_ext},
    scripts=['scripts/gat-run.py',
             'scripts/gat-great.py',
             'scripts/gat-compare.py',
             'scripts/gat-plot.py'],
    zip_safe=False,
)
