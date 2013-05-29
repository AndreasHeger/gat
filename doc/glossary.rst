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
      
   sample
      one set of :term:`sampled segments`.

   isochore
      Usually, isochores are defined as genomic regions of homogeneous
      G+C content. In GAT, isochores can mean any genomic regions with
      a shared property. Isochores are used to eliminate the effect
      of confounders in an analysis.

   isochores
      plural of :term:`isochore`.

   isochore workspaces
      the workspace split according to :term:`isochores`.

   bed
      an interval format. Intervals are in a tab-separated list denoted
      as ``contig``, ``start``, ``end``. A bed file can contain several
      :term:`tracks` which will be treated independently. Tracks are either delineated
      by a line of the format: ``track name=<trackname>`` or by an optional 
      fourth column (field ``name``). See UCSC_ for more information about
      bed files.

   fold 
      fold change, computed as the ratio of observed over expected

   p-value
      significance of association. The P-value is an estimate of the
      probability to obtain an observed (or larger) overlap between 
      two segment sets by chance.

   q-value
      multiple testing corrected p-value. The qvalue can either be
      computed using the q-value procedure suggested by 
      `Storey et al. (2002)`_ or using other FDR approaches
      such as Bonferroni, Benjamini-Hochberg, etc. implemented
      in the R_ function ``p.adjust``.

   pvalue
      the column containing the :term:`p-value`.

   qvalue
      the column containing the :term:`q-value`.

   track
      a set of segments within a :term:`bed` formatted file. Tracks
      are either grouped using the ``track name=<trackname>`` prefix
      or by an optional fourth column (the ``name`` column).

   annotation
      a genomic annotation, which is a collection of intervals.

   observed
      the observed value of a metric.

   expected
      the expected value of a metric as determined by simulations.

   fold
      the fold change, defined as the ratio of observed over expected.

   l2fold
      the logarithm (base 2) of :term:`fold` change.

   tracks
      plural of :term:`track`

   samples
      a set of samples (see :term:`sample`)

.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _R: http://www.r-project.org
.. _Storey et al. (2002): http://genomics.princeton.edu/storeylab/papers/directfdr.pdf
