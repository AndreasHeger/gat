==========
Background
==========

History
=======

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

.. Use of GAT in studies

.. GAT has been used in the following published studies:

.. _Ponjavic et al (2007): http://genome.cshlp.org/content/17/5/556.short
.. .. _Heger at al. 

Comparison to other methods
===========================

Testing for the association between genomic features is a field of
long-standing interest in genomics and has gained considerable
traction with the publication of large scale genomic data sets such as
the ENCODE_ data.

Generally we believe that the problem of testing for association has
not been fully resolved and advise every genomicist to apply several
methods. The list of tools/services below is not exhaustive:

GREAT_ (`MacLean et al. (2010)`_) uses a binomial test to test if
transcription factor binding sites are associated with regulatory
domains of genes. GREAT_ has a convenient web interface with many
annotations pre-loaded. Compared to GREAT_, GAT can measure depletion
and can 

GenometriCorr_ (`Favorov et al. (2012)`_) compute a variety of distance
metrics when comparing two interval sets and then apply a set of
statistical tests to measure association. GenomicCorr is a good
exploratory tool to generate hypotheses about the relationships
of two genomic sets of intervals. Compared to GenometriCorr_, GAT
can simulate more realistic genomic scenarios, for example, segments might
not occur in certain regions (due to mapping problems) or occur at
reduced frequency (G+C biases). 

The GSC_ (`The Encode Project Consortium (2012)`_) metric (for Genome
Structure Correlation) is inspired by the analysis of approximately piecewise
stationary time series . The GSC metric estimates the significance of
an association metric by estimating the random expectation of the
association metric using randomly chosen intervals on the genome. This
expectation is then used to test if the observed value of the metric
(nucleotide overlap, region overlap, ...) is higher than expected. 
The method is described in the supplemental details of the first 
and recent ENCODE_ papers (`Birney et al. (2007)_`, ...) and 
`here <http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.aoas/1294167794>`_.


BITS_ (`Layer et al. (2013)`_) (Binary Interval Search) is a method
to perform quick overlap queries between genomic data sets. It
implements a Monte-Carlo method for simulation that is particularly
suited towards making all on all comparisons between a large number 
data sets.

.. _GREAT: http://bejerano.stanford.edu/great/public/html/
.. _MacLean et al. (2010): http://www.ncbi.nlm.nih.gov/pubmed/20436461
.. _GenometriCorr: http://genometricorr.sourceforge.net/
.. _Favorov et al. (2012): http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002529#pcbi-1002529-g001
.. _GSC: http://www.encodestatistics.org/
.. _The Encode Project Consortium (2012): http://www.nature.com/nature/journal/v489/n7414/full/nature11247.html
.. _BITS: https://github.com/arq5x/bits
.. _Layer et al. (2013): http://www.ncbi.nlm.nih.gov/pubmed/23129298

Benchmark
---------

We used the example from :ref:`Tutorial1` to perform a rough
comparison between various methods. In all cases, we used
``n = 1000`` for simulations. Times are wall-clock times.
Please note that this is not a rigorous benchmark.

+------+------+------------+--------+-----------+--------+-----+
|Method|Set1  |Set2        |Observed|Expected   |P-value |Time |
+------+------+------------+--------+-----------+--------+-----+
|BITS_ |srf   |jurkat      |450     |5.24       |<0.001  |43s  |
+------+------+------------+--------+-----------+--------+-----+
|BITS_ |srf   |hepg2       |381     |9.87       |<0.001  |39s  |
+------+------+------------+--------+-----------+--------+-----+
|BITS_ |srf   |hepg2/jurkat|9       |5.7        |0.13    |28s  |
+------+------+------------+--------+-----------+--------+-----+
|BITS_ |jurkat|hepg2       |47237   |3548       |<0.001  |106s |
+------+------+------------+--------+-----------+--------+-----+
|GSC_  |srf   |jurkat      |        |           |0.0004  |58s  |
+------+------+------------+--------+-----------+--------+-----+
|GSC_  |srf   |hepg2       |        |           |6.9E-11 |54s  |
+------+------+------------+--------+-----------+--------+-----+
|GSC_  |srf   |hepg2/jurkat|        |           |5.23E-7 |40s  |
+------+------+------------+--------+-----------+--------+-----+
|GSC_  |jurkat|hepg2       |        |           |0       |159s |
+------+------+------------+--------+-----------+--------+-----+
|GAT   |srf   |jurkat      |20183   |247.6      |<0.001  |11s  |
+------+------+------------+--------+-----------+--------+-----+
|GAT   |srf   |hepg2       |18965   |601.4      |<0.001  |11s  |
+------+------+------------+--------+-----------+--------+-----+
|GAT   |srf   |hepg2/jurkat|425     |327.3      |0.21    |11s  |
+------+------+------------+--------+-----------+--------+-----+
|GAT   |jurkat|hepg2       |6163503 |457332.8   |<0.001  |316s |
+------+------+------------+--------+-----------+--------+-----+

BITS_ and GAT are fairly comparable, even though they use different
metrics for the association (number of segments overlapping versus
number of nucleotides overlapping). GAT is quicker on smaller data
sets, while BITS_ outperforms on large datasets.

GSC_ reports a significant association in the comparison between
srf and dhs intervals specific to hepg2 cells, while the other two
tools do not, which is the biologically plausible result. It is
difficult to say if there indeed is an association, or GSC_ is 
overestimating association.

.. _ENCODE: http://genome.ucsc.edu/ENCODE/
.. _Birney et al. (2007): http://www.nature.com/nature/journal/v447/n7146/extref/nature05874-s1.pdf
