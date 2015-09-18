The genomic association tester (GAT) tests for association between
intervals in a genome analysis context.

GAT answers the question, if a set of intervals on a genome is
statistically significantly associated with one or more other sets of
genomic intervals. Assocation typically means overlap, but also other
measures like the distance between elements can be tested.

The tests are performed by simulation and within a genomic
context. The simulation part takes care of segment length
distribution. The genomic context takes into account chromosomal
location and optionally isochores.

The software is available for download on pypi:

    pip install gat

Please see http://gat.readthedocs.org/en/latest/ for a manual and tutorials.

A GAT user group is at:

https://groups.google.com/forum/?fromgroups#!forum/gat-user-group

The latest release can be downloaded from pypi at
https://pypi.python.org/pypi/gat
