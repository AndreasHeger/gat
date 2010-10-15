========
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
