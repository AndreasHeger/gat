.. gat documentation master file, created by
   sphinx-quickstart on Mon Jul 12 11:11:43 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================
Genomic Association Tester (GAT)
================================

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

Usage instructions
==================

The *gat* tool is controlled via the :file:`gatrun.py` script. This 
script requires the following input:

   1. A set of :term:`intervals` ``S`` with :term:`segments` to test.
   2. A set of :term:`intervals` ``A`` with :term:`annotations` to test against.
   3. A set of :term:`intervals` ``W`` describing a :term:`workspace` 

The input can be in a variety of formats, for example :term:`bed`
formatted files. The script is then run as::

   gatrun.py 
      --segment-file=segments.bed.gz 
      --workspace-file=workspace.bed.gz 
      --annotation-file=annotations_architecture.bed.gz  

The script recognizes *gzip* compressed files by the suffix ``.gz``.

The principal output is a tab-separated table of pairwise comparisons between
each segment and track. The table will be written to stdout, unless the option 
``--stdout`` is redirected to a file. The columns in the table are:

   track 
      the segment track
   annotation
      the annotation track
   observed        
      the observed count
   expected        
      the expected count based on the samples
   CI95low 
      the value at the 5% percentile of sample
   CI95high        
      the value at the 95% percentile of sample
   stddev  
      the standard deviation of samples
   fold    
      the fold enrichment, given by the ratio observed / expected
   pvalue
      the pvalue of enrichment/depletion

Further output files go to files named according to ``--filename-output-pattern``.
``%s`` are substituted with section names.

*Count* here denotes the measure of association and defaults
to *number of overlapping nucleotides*.

Submitting multiple files
-------------------------

All of the options *--segment-file*, *--workspace-file*, *--annotation-file* 
can be used several times on the command line. What happens with multiple files
depends on the file type:

   1. Multiple *--segment-file* entries are added to the list of segment lists
      to test with.

   2. Multiple *--annotation-file* entries are added to the list of annotation lists
      to test against.

   3. Multiple *--workspace* entries are intersected to create a simple workspace.

Generally, *gat* will test *m* segment lists against *n* annotation lists in all
*m * n* combinations.

Adding isochores
----------------
   
Isochores are genomic segments with common properties that are potentially correlated
with the segments of interest and the annotations, but the correlation is not of interest
here. For example, consider a CHiP-Seq experiment and the testing if CHiP-Seq intervals
are close to genes. G+C rich regions in the genome are gene rich, while at the same time 
there is possibly a nucleotide composition bias in the CHiP-Seq protocol depleting A+T
rich sequence. An association between genes and CHiP-Seq intervals might simply be due
to the G+C effect. Using isochores can control for this effect to some extent.

Isochores split the :term:`workspace` into smaller workspaces of similar properties,
so called *isochore workspaces*. Simulations are performed for each :term:`isochore workspaces` 
separately. At the end, results for each all isochore workspaces are aggregated.

In order to add isochores, us the *--isochore-file* command line option.

Choosing measures of association
--------------------------------

Counters describe the measure of association that is tested. Counters
are selected with the command line option ``--counter``. Available 
counters are:

   1. ``nucleotide-overlap``: number of bases overlapping
   2. ``segment-overlap``: number of intervals intersecting

Changing the PValue method
--------------------------

By default, *gat* returns the empirical pvalue based on the sampling
procedure. The minimum :term:`pvalue` is ``1 / number of samples``.

If the option ``--pvalue`` is set to ``--pvalue=norm``, pvalues are
computed by fitting a normal distribution to the samples.

Multiple testing
----------------

*gat* can use the procedure by `Storey et al. (2002)`_ to compute a
:term:`qvalue` for each pairwise comparison. The implementation
is equivalent to the qvalue_ package implemented in R_.

Caching
-------

*gat* can save and retrieve samples from a cache ``--cache=cache_filename``.
If :file:`cache_filename` does not exist, samples will be saved to the
cache after computation. If :file:`cache_filename` does already exist,
samples will be retrieved from the cache instead of being re-computed.
Using cached samples is useful when trying different :term:`Counters`.

If the option ``--counts-file`` is given, *gat* will skip the sampling
and counting step completely and read observed counts from 
``--count-file=counts_filename``.

List of all command-line options
================================



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


Glossary
========

.. glossary::
  
   Interval
      a (genomic) segment with a chromosomal location (contig,start and end).

   Intervals
      a set of one or more intervals.

   workspace
      the genomic regions accessible for simulation.

   annotations
      sets of intervals annotating various regions of the genome.

   segments
      sets of intervals whose association is tested with :term:`annotations`.

   bed
      an interval format. Intervals are in a tab-separated list denoted
      as ``contig``, ``start``, ``end``. A bed file can contain several
      tracks which will be treated independently. Tracks are either delineated
      by a line of the format: ``track name=<trackname>`` or by an optional 
      fourth column (field ``name``). See UCSC_ for more information about
      bed files.

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

Background
==========

This module has been inspired by the TheAnnotator tool first
used in the analysis by `Ponjavic et al (2007)`_ by Gerton Lunter and early 
work of Caleb Webber.

The differences are:

   * permit more measures of association. The original Annotator used nucleotide 
     overlap, but other measures might be useful like number of elements overlapping
     by at least x nucleotides, proximity to closest element, etc.
 
   * easier user interface and using standard formats.

   * permit incremental runs. Annotations can be added without recomputing the samples.

   * faster.


.. _qvalue: http://genomics.princeton.edu/storeylab/qvalue/linux.html
.. _R: http://www.r-project.org
.. _Ponjavic et al (2007): http://genome.cshlp.org/content/17/5/556.short
.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _Storey et al. (2002): http://genomics.princeton.edu/storeylab/papers/directfdr.pdf
