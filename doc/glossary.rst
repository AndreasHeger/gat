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

   segments of interest
      sets of intervals whose association is tested with
      :term:`annotations`

   sampled segments
      randomized versions of the :term:`segments of interest`. A set
      of sampled segments is one :term:`sample`.
      
   sample:
      on set of :term:`sampled segments`.

   bed
      an interval format. Intervals are in a tab-separated list denoted
      as ``contig``, ``start``, ``end``. A bed file can contain several
      :term:`tracks` which will be treated independently. Tracks are either delineated
      by a line of the format: ``track name=<trackname>`` or by an optional 
      fourth column (field ``name``). See UCSC_ for more information about
      bed files.

   fold 
      fold change, computed as the ratio of observed over expected

   pvalue
      significance of association. The P-value is an estimate of the
      probability to obtain an observed (or larger) overlap between 
      two segment sets by chance.

   track
      a set of segments within a :term:`bed` formatted file. Tracks
      are either grouped using the ``track name=<trackname>`` prefix
      or by an optional fourth column (the ``name`` column).

   samples
      a set of samples (see :term:`sample`)

.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat#format1
