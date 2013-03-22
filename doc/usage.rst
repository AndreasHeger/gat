.. _Usage:

================================
Usage instructions
================================

This page describes basic and advanced usage of GAT.

A list of all command-line options is available via::
 
   gat-run.py --help

Basic usage
===========

The *gat* tool is controlled via the :file:`gat-run.py` script. This 
script requires the following input:

   1. A set of :term:`intervals` ``S`` with :term:`segments of interest` to test.
   2. A set of :term:`intervals` ``A`` with :term:`annotations` to test against.
   3. A set of :term:`intervals` ``W`` describing a :term:`workspace` 

GAT requires :term:`bed` formatted files. In its simplest form, GAT is then run as::

   gat-run.py 
      --segment-file=segments.bed.gz 
      --workspace-file=workspace.bed.gz 
      --annotation-file=annotations.bed.gz  

The script recognizes *gzip* compressed files by the suffix ``.gz``.

The principal output is a tab-separated table of pairwise comparisons between
each :term:`segments of interest` and :term:`annotations`. The table
will be written to stdout, unless the option ``--stdout`` is given
with a filename to which output should be redirected.

The main columns in the table are:

   track 
      the :term:`segments of interest` :term:`track`
   annotation
      the :term:`annotations` :term:`track`
   observed        
      the observed count
   expected        
      the expected count based on the :term:`sampled segments`
   CI95low 
      the value at the 5% percentile of :term:`samples`
   CI95high        
      the value at the 95% percentile of :term:`samples`
   stddev  
      the standard deviation of :term:`samples`
   fold    
      the fold enrichment, given by the ratio observed / expected
   l2fold
      log2 of the fold enrichment value
   pvalue_
      the pvalue_ of enrichment/depletion
   qvalue
      a multiple-testing corrected :term:`pvalue`. See 
      `multiple testing correction`_.

Additionally, there are the following columns:
    track_nsegments 
       number of segments in :term:`track` in :term:`segments of interest` 
    track_size 
       number of residues in covered by :term:`track` in :term:`segments of interest` within the :term:`workspace`
    track_density
       fraction of residues in :term:`track` in :term:`segments of interest` within the :term:`workspace`
    annotation_nsegments
       number of segments in :term:`track` in :term:`annotations`.
    annotation_size
       number of residues in covered by :term:`track` in :term:`annotations` within the :term:`workspace`
    annotation_density
       number of residues in covered by :term:`track` in :term:`annotations` within the :term:`workspace`
    overlap_nsegments       
       number of segments in overlapping between :term:`segments of interest` and :term:`annotations`
    overlap_size    
       number of nucleotides overlapping between :term:`segments of interest` and :term:`annotations`
    overlap_density
       fraction of residues overlapping between :term:`segments of interest` and :term:`annotations`
       within :term:`workspace`
    percent_overlap_nsegments_track 
       percentage of segments in :term:`segments of interest` overlapping :term:`annotations`
    percent_overlap_size_track
       percentage of nucleotides in :term:`segments of interest` overlapping :term:`annotations`
    percent_overlap_nsegments_annotation
       percentage of segments in :term:`annotations` overlapping :term:`segments of interest`
    percent_overlap_size_annotation 
       percentage of nucleotides in :term:`annotations` overlapping :term:`segments of interest`
    description
       additional description of track (requires ``--descriptions`` to
       be set).

Further output files such as auxiliary summary statistics go to files
named according to ``--filename-output-pattern``. The argument to
``filename-output-pattern`` should contain one ``%s`` placeholder, 
which is then substituted with section names.

*Count* here denotes the measure of association and defaults to *number of overlapping nucleotides*.

Advanced Usage
==============

Submitting multiple files
-------------------------

All of the options *--segment-file*, *--workspace-file*, *--annotation-file* 
can be used several times on the command line. What happens with multiple files
depends on the file type:

   1. Multiple *--segment-file* entries are added to the list of
      :term:`segments of interest` to test with.

   2. Multiple *--annotation-file* entries are added to the list of :term:`annotations`
      to test against. 

   3. Multiple *--workspace* entries are intersected to create a single workspace.

Generally, *gat* will test *m* :term:`segments of interest` lists
against *n* :term:`annotations` lists in all *m * n* combinations.

Within a :term:`bed` formatted file, different :term:`tracks` can
be separated using a UCSC formatted ``track`` line, such as this::

   track name="segmentset1"
   chr1	 23	100
   chr3	 50	2000
   track name="segmentset2"
   chr1	 1000	2000
   chr3	 4000	5000

or alternatively, using the fourth column in a :term:`bed` formatted file::
	
   chr1	23	100	segmentset1
   chr3	50	2000	segmentset1
   chr1	1000	2000	segmentset2
   chr3	4000	5000	segmentset2

The latter takes precedence. The option `--ignore-segment-tracks``
forces gat to ignore the fourth column and consider all intervals
to be from a single interval set.

.. note::
   Be careful with bed-files where each interval gets a unique
   identifier. Gat will interprete each interval as a separate
   segment set to read. This is usually not intended and causes
   gat to require a very large amount of memory.
   (see the option ``--ignore-segment-tracks``

By default, tracks can not be split over multiple files. The option
``--enable-split-tracks`` permits this.

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

In order to add isochores, use the *--isochore-file* command line option.

Choosing measures of association
--------------------------------

Counters describe the measure of association that is tested. Counters
are selected with the command line option ``--counter``. Available 
counters are:

1. ``nucleotide-overlap``: number of bases overlapping [default]

2. ``segment-overlap``: number of intervals intervals in the
   :term:`segments of interest` overlapping :term:`annotations`. A single
   base-pair overlap is sufficient.

3. ``segment-mid-overlap``: number of intervals in the
   :term:`segments of interest` overlapping at their midpoint
   :term:`annotations`.

4. ``annotations-overlap``: number of intervals in the
   :term:`annotations` overlapping :term:`segments of intereset`. A single
   base-pair overlap is sufficient.

5. ``segment-mid-overlap``: number of intervals in the
   :term:`annotations` overlapping at their midpoint 
   :term:`segments of intereset`
   
Multiple counters can be given. If only one counter is provided, the
output will be to stdout. Otherwise, separate output files will be
created each counter. The filename can be controlled with the
``--output-table-pattern`` option.

Changing the PValue method
--------------------------

By default, *gat* returns the empirical :term:`pvalue` based on the sampling
procedure. The minimum :term:`pvalue` is ``1 / number of samples``.

Sometimes the lower bound on p-values causes methods that estimate the
FDR to fail as the distribution of p-values is atypical. In order to
estimate lower pvalues, the number of samples needs to be increased.
Unfortunately, the run-time of gat is directly proportional to the number of
samples.

A solution is to set the option ``--pvalue-method`` to ``--pvalue-method=norm``. In that
case, pvalues are estimated by fitting a normal distribution to the
samples. Small p-values are obtained by extrapolating from this fit.

Multiple testing correction
---------------------------

*gat* provides several methods for controlling the `false discovery
rate`_. The default is to use the Benjamini-Hochberg procedure.
Different methods can be chosen with the ``--qvalue-method`` option.

``--qvalue-method=storey`` uses the procedure by `Storey et al. (2002)`_ to compute a
:term:`qvalue` for each pairwise comparison. The implementation
is in its functionality equivalent to the qvalue_ package implemented 
in R_.

Other options are equivalent to the methods as implemented in the
R_ function ``p.adjust``.

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

Using multiple CPU/cores
========================

GAT can make use of several available CPU/cores if available. Use
the ``--num-threads=#`` option in order to specify how many CPU/cores
GAT will make use of. The default ``--num-threads=0`` means that GAT
will not use any multiprocessing.

Outputting intermediate results
-------------------------------

A variety of options govern the output of intermediate results by gat.
These options usually accept patterns that represent filenames with 
a ``%s`` as a wild card character. The wild card is replaced with
various keys. Note that the amount of data output can be substantial.

``--output-counts-pattern``
   output counts. One file is created for each counter. Counts output files
   are required for :ref:`gat-compare`.   

``--output-plots-pattern``
   create plots (requires matplotlib_). One plot for each annotation
   is created showing the distribution of expected counts and the
   observed count. Also, outputs the distribution of p-values and
   q-values.
   
``--output-samples-pattern``
   output :term:`bed` formatted files with individual samples.

Other tools
===========

.. _gat-compare:

gat-compare
-----------

The gat-compare tool can be used to test if the fold changes found in
two or more different gat experiments are significantly different from
each other.

This tool requires the output files with counts created using the
``--output-counts-pattern`` option.

For example, to compare if fold changes are signficantly different
between two cell lines, execute::

   gat-run.py --segments=CD4.bed.gz <...>
   --output-counts-pattern=CD4.%s.overlap.counts.tsv.gz
   gat-run.py --segments=CD14.bed.gz <...>
   --output-counts-pattern=CD14.%s.overlap.counts.tsv.gz

   gat-compare.py CD4.nucleotide-overlap.counts.tsv.gz CD14.nucleotide-overlap.counts.tsv.gz

.. _gat-plot:

gat-plot
--------

Plot gat results.

.. _gat-great:

gat-great
---------

Perform a GREAT_ analysis::

   gat-great.py 
      --segment-file=segments.bed.gz 
      --workspace-file=workspace.bed.gz 
      --annotation-file=annotations.bed.gz  

.. _pvalue: http://en.wikipedia.org/wiki/P-value
.. _qvalue: http://genomics.princeton.edu/storeylab/qvalue/linux.html
.. _R: http://www.r-project.org
.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _Storey et al. (2002): http://genomics.princeton.edu/storeylab/papers/directfdr.pdf
.. _false discovery rate: http://en.wikipedia.org/wiki/False_discovery_rate
.. _matplotlib: http://matplotlib.org/
.. _GREAT: http://bejerano.stanford.edu/great/public/html/
