================================
Genomic Association Tester (GAT)
================================

Welcome to the home page of the Genomic Association Tester (*GAT*).

Overview
========

A common question in genomic analysis is whether two sets of genomic
intervals overlap significantly. This question arises, for example, in
the interpretation of ChIP-Seq or RNA-Seq data. Because of complex
genome organization, its answer is non-trivial.

The Genomic Association Tester (GAT) is a tool for computing
the significance of overlap between multiple sets of genomic
intervals. GAT estimates significance based on simulation and
can take into account genome organization like isochores and correct for
regions of low mapability.

GAT accepts as input standard genomic file formats and can be used in
large scale analyses, comparing the association of multiple sets of
genomic intervals simultaneously. Multiple testing is controlled using
the false discovery rate.

In this manual, the :ref:`Introduction` covers the basic concepts of
GAT. In order to get an idea of typical use cases, see the
:ref:`Tutorials` section. The :ref:`Usage` section contains
a complete usage reference.

Contents
========

.. toctree::
   :maxdepth: 2

   introduction.rst   
   installation.rst
   tutorials.rst
   usage.rst
   fold.rst
   performance.rst
   background.rst
   glossary.rst
   release.rst
   
Developers' notes
====================

The following section contains notes for developers.

.. toctree::
   :maxdepth: 2

   notes.rst   
   testing.rst
   simulators.rst

..   examples.rst
..   implementation.rst
..   tutorial.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. .. todo::

..    * permit incremental runs
..    * caching of segment size distribution
..    * pre-compute workspace/segment intersection
..    * parallel computation of samples 
..    * profile and tune
..    * add GFF/GTF support
..    * add additional counters:
..       Distance between segment and annotation
..       1. Closest distance of segment to annotation
..       2. Closest distance of annotation to segment
..       3. Fix Segment density counter
..    * benchmarking - add diagrams for experiment setup 
..    * benchmarking - number of segment distributions
..    * compute summarry statistics on workspaces, segments and
..      annotations in order to detect excessive fragmentation.
..    * Benchmark this.

..    *  Summary parameters
..    *  Parameters are: 
..        * number of segments, number of nucleotides
..        * number of segments split by a boundary, nucleotides extending outside boundary
..        * segment size distribution untruncated (min, max, mean, media, q1, q3)
..        * segment size distributino truncated

..      Count segments and annotations
..      Count per contig in the workspace, count per isochore

..      Optionally: do all of this for sampled segments


..    * add galaxy interface.
..    * reorganize package structure
..    * check pip installation
..    * check requirements   

  
