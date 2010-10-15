.. gat documentation master file, created by
   sphinx-quickstart on Mon Jul 12 11:11:43 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================
Usage instructions
================================

Basic usage
-----------

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

Multiple testing correction
---------------------------

*gat* can use the procedure by `Storey et al. (2002)`_ to compute a
:term:`qvalue` for each pairwise comparison. The implementation
is in its functionality equivalent to the qvalue_ package implemented 
in R_.

Caching sampling results
------------------------

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

A list of all command-line options is available via::
 
   gatrun.py --help

.. _qvalue: http://genomics.princeton.edu/storeylab/qvalue/linux.html
.. _R: http://www.r-project.org
.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _Storey et al. (2002): http://genomics.princeton.edu/storeylab/papers/directfdr.pdf
