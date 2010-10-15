===============
Implemenation
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




