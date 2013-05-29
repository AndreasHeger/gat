===============
Implementation
===============

Caching samples
===============

Several options are able to cache samples:

   * using pytables: permits compression and easy access
       to results through its hierarchic structure that maps
       straight onto a nested dictionary.
       Slowish I/O as SegmentList does not map directly to pytables
       types and values need to be committed individually. Easy
       to code, permits extending.

   * compressed :term:`bed` files. They are a natural format in
       genomic context and would permit direct visualization in
       genome browsers. Slightly wasteful because of contigs
       being names.

   * binary segment lists. Fast I/O, but more coding required.
       Requires separate management of tracks, samples and 
       contigs.

To illustrate the space usage, I consider an experiment with
17 tracks and 59193 intervals. The compressed input ".bed.gz"
data takes 496k total. For 10,000 samples, the size of the cached
samples would be:
   
   * pytables: 3h and 11Gb (experimentally determined)
   * binary block compressed segment list: 1h15' and 3.3Gb (experimentally determined,
      without index). Reading the data instead of sampling takes 10'. This includes	
      the time for counting.
   * compressed :term:`bed` files: 10,000 * 496k = 4.96Gb.
   * binary uncompressed segment lists: 10,000 * 59,193 * 8bytes = 4.7Gb
     (without index)

Without caching, sampling and p-value computation takes 1h2'.

Manual coding is still the best.

GREAT
======

The GREAT algorithm works as follows within the GAT
terminology. In GREAT, annotations are overlapping regulatory domains
of genes with GO terms againts which segments of interest are tested.

For each annotation A:

1. Compute total number of bases annotated with annotation A: b_A
2. Compute total number of bases in workspace: b_W
3. Compute fraction of 1 and 2: p_A = fraction of genome annotated
   with A = b_A / b_W
2. Compute number of times a segment hits a certain annotation: k_A
3. Compute number of times a segment is in the workspace: n_S

P = Pr_binom( k >= k_A | n = n_S, p = p_A )

Testing for fold change differences
===================================

Testing for fold change differences using the fold change ratios of
samples and tests them against 1.

It is important to perform the test on the fold changes and not the distributions
of counts. Counts between different combinations of segments and
annotations are not comparable. A larger set of :term:`segments of
interest` will increase the expected counts compared to a smaller set.
Similarly, a set with larger segments will have a higher expected
overlap than a smaller set. Thus, only after normalization with the
observed value the comparison can be made.

I use the ratio of observed values, but I believe the difference could be used
as well. Using the ratio in log-space facilitates plotting. It might
also better compensate for extreme values if the expected overlap is
small and there are many samples with no overlap.











