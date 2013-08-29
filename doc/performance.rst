.. _performance:

===========
Performance
===========

Memory usage
============

GAT is limited by the amount of memory available. Factors affecting
GAT's memory usage are:

* the number of segments. Memory requirements grow linearly with
  the number of sets and the number of intervals in the sets
  :term:`segments of interest`, :term:`annotations` and
  :term:`workspace`. The total number of segments 

* the number of samples. For each sample, statistics on nucleotide
  overlap are being kept. Thus memory requirements grow linearly with
  the number of samples.

If memory consumption is problematic, the following might help:

* Only working with one set of :term:`segments of interest`.
* Only working with few sets of :term:`annotations`.
* In a pairwise comparison, use the smaller set as :term:`segments of
  interest` and the larger as :term:`annotations`.
* Avoiding fragmentation of the workspace.
* Caching of samples

Memory consumption is usually not problematic and for typical sets
reaches a few Gigabytes. If memory consumption explodes it is usually
a problem with the input adat. 

GAT has not been optimized for memory usage. If memory consumption is
a problem, please contact the developers.

Runtime performance
===================

As with memory usage, run-time performance of GAT is linearly related
to the number of segments and the number of samples.

* the number of segments in the :term:`segments of interest`. More
  segments require more sampling steps. Usually, an input size with
  twice as many segments will take twice as long. Runtime behaviour might be
  worse in extreme cases where the sampler has difficulties placing the last segments, for
  example if the segment density is high.

* the number of samples. Each additional sample will require an
  additional amount of time.

* the number of segments in the :term:`annotations`. Usually less of a
  factor, but it becomes a factor if overlap statisticts need to be
  computed for many different annotations. In this case, generating a
  single sample might be quicker than computing overlap statisticts of
  that sample with 

As runtime performance is a linear function of most variables, runtime
can be reduced by using multiple CPUs or cores (see
:ref:`multiplecores`).

Sampling performance
====================

In order to check for biases in the sampling procedure, we run
automated tests to check for even coverage of nucleotides by the
sampler and absence of edge effects.

.. toctree::
   :maxdepth: 2

   testingSamplerSamplerAnnotator.rst

