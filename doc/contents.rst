================================
Genomic Association Tester (GAT)
================================

Welcome to the home page of the Genomic Association Tester (*GAT*).

Overview
========

The genomic association tester (GAT) tests for association
between sets of genomic intervals.

GAT tests if two sets of genomic intervals are associated more than expected by chance.
Assocation typically means nucleotide overlap, but other 
measures such as the distance between elements or the number of
overlapping segments can be used. 

The tests are performed by simulation within a genomic context. 
The simulation part takes care of segment length distribution. 
The genomic context takes into account chromosomal location and optionally 
isochores.

Contents
========

.. toctree::
   :maxdepth: 2

   tutorial.rst
   usage.rst
   examples.rst
   background.rst
   testing.rst
   implemenation.rst
   glossary.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. todo::

   * permit incremental runs
   * caching of segment size distribution
   * pre-compute workspace/segment intersection
   * parallel computation of samples 
   * profile and tune
   * add GFF/GTF support
   * add additional samplers
   * add additional counters:
      Distance between segment and annotation
      1. Closest distance of segment to annotation
      2. Closest distance of annotation to segment
      3. Fix Segment density counter

