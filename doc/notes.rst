=====
Notes
=====

Sampling strategies
===================

Sampling creates a new set of :term:`interval` ``P``. There
are several different strategies possible.

Annotator strategy
------------------

In the original Annotator strategy, samples are created in a two step procedure. 
First, the length of a sample segment is chosen randomly from the empirical segment 
length distribution. Then, a random coordinate is chosen. If the
sampled segment does not overlap with the workspace it is rejected. 

Before adding the segment to the sampled set, it is truncated to 
the workspace.

If it overlaps with a previously sampled segment, the segments
are merged. Thus, bases shared between two segments are not counted 
twice.

Sampling continues, until exactly the same number of bases overlap between
the ``P`` and the ``W`` as do between ``S`` and ``W``.
      	 
Note that the length distribution of the intervals in ``P`` might be different 
from ``S``.

The length is sampled from the empirical histogram of segment sizes. The
granularity of the histogram can be controlled by the options ``-bucket-size``
and ``--nbuckets``. The largest segment should be smaller than ``bucket-size * nbuckets``.
Increase either if you have large segments in your data set, but smaller
values of ``nbuckets`` are quicker.

This method is quick if the workspace is large. If it is small, 
a large number of samples will be rejected and the procedure 
slows down.

This sampling method is compatible with both distance and overlap
based measures of associaton. 

Workspaces and isochores
++++++++++++++++++++++++

Workspaces limit the genomic space available for segments and annotations.
Isochores split a workspace into smaller parts that permit to control for
confounding variables like G+C content.

The simplest workspace is the full genome. For some analyses it might be better 
to limit to analysis to autosomes. 

Examples for the use of isochores could be to analyze chromosomes or chromosomal arms
separately. 

If isochores are used, the length distribution and nucleotide overlaps are counted per isochore
ensuring that the same number of nucleotides overlap each isochore in ``P`` and ``S`` and the
length distributions per isochore are comparable. 

Empirical length distribution
+++++++++++++++++++++++++++++

The empirical length distribution is created from all :term:`intervals`
in ``S``. The full segment length is chosen even if there is partial overlap.
Optionally, the segment can be truncated. From Gerton::

   What is the best choice depends on the data. Not truncating can lead 
   to a biased length distribution if it is expected that segments that 
   only partially overlap the workspace have very different lengths. However, 
   truncating can lead to spurious short segments.
