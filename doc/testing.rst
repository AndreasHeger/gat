============
Benchmarking
============

This section contains quality control plots from the unit testing.

Sampling
========

The following plots benchmark the segment sampling behaviour of GAT and Annotator.
The tests below sample from a segment list *10* segments with each segment of size *100*.

.. toctree::
   :maxdepth: 2

   testingSamplerGatSegments.rst
   testingSamplerGatAnnotator.rst
   testingSamplerGatUniform.rst
   testingSamplerAnnotator.rst
   testingPosition.rst

Statistics
==========

For 1-sized fragments (i.e. SNPs), the statistics can be checked against
a hypergeometric distribution (sampling without replacement). All the
tests below use a single continuous workspace of 1000 nucleotides seeded
with a varying number of SNPs. 

.. figure:: ../test/test_testSingleSNP.TestStatsSNPSampling.png
   :width: 500

   Test with a single SNP. Here, there are no issues with multiple hits.
   The workspace contains a single annotation of increasing size (1,3,5,...,99)

.. figure:: ../test/test_testMultipleSNPsFullOverlap.TestStatsSNPSampling.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains a
   single annotation of size (10, 15, ..., 105). All SNPs overlap the annotated
   part of the workspace and hence all results are highly signficant. 

.. figure:: ../test/test_testMultipleSNPsPartialOverlap.TestStatsSNPSampling.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains 
   a single annotation. Annotations are all of size *10*, but the overlap of
   SNPs with annotations varies from 0 to 10.

Statistics
==========

Gat
---

SNPs
++++

For 1-sized fragments (i.e. SNPs), the statistics can be checked against
a hypergeometric distribution (sampling without replacement). All the
tests below use a single continuous workspace of 1000 nucleotides seeded
with a varying number of SNPs. 

.. figure:: ../test/testSingleSNP.TestStatsGat.png
   :width: 500

   Test with a single SNP. Here, there are no issues with multiple hits.
   The workspace contains a single annotation of increasing size (1,3,5,...,99)

.. figure:: ../test/testMultipleSNPsFullOverlap.TestStatsGat.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains a
   single annotation of size (10, 15, ..., 105). All SNPs overlap the annotated
   part of the workspace and hence all results are highly signficant. 

.. figure:: ../test/testMultipleSNPsPartialOverlap.TestStatsGat.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains 
   a single annotation. Annotations are all of size *10*, but the overlap of
   SNPs with annotations varies from 0 to 10.

.. figure:: ../test/testWorkspaces.TestStatsGat.png
   :width: 500

   workspace = 500 segments of size 1000, separated by a gap of 1000
   annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
   segments = a SNP every 100 bp

Intervals
+++++++++

.. figure:: ../test/testIntervalsPartialOverlap.TestStatsGat.png
   :width: 500

   In this test, there is one segment of size *100*. Annotations
   are of size *100* with decreasing overlap.


Annotator
---------

SNPs
++++

For 1-sized fragments (i.e. SNPs), the statistics can be checked against
a hypergeometric distribution (sampling without replacement). All the
tests below use a single continuous workspace of 1000 nucleotides seeded
with a varying number of SNPs. 

.. figure:: ../test/testSingleSNP.TestStatsTheAnnotator.png
   :width: 500

   Test with a single SNP. Here, there are no issues with multiple hits.
   The workspace contains a single annotation of increasing size (1,3,5,...,99)

.. figure:: ../test/testMultipleSNPsFullOverlap.TestStatsTheAnnotator.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains a
   single annotation of size (10, 15, ..., 105). All SNPs overlap the annotated
   part of the workspace and hence all results are highly signficant. 

.. figure:: ../test/testMultipleSNPsPartialOverlap.TestStatsTheAnnotator.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains 
   a single annotation. Annotations are all of size *10*, but the overlap of
   SNPs with annotations varies from 0 to 10.

.. figure:: ../test/testWorkspaces.TestStatsTheAnnotator.png
   :width: 500

   workspace = 500 segments of size 1000, separated by a gap of 1000
   annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
   segments = a SNP every 100 bp

Intervals
+++++++++

.. figure:: ../test/testIntervalsPartialOverlap.TestStatsTheAnnotator.png
   :width: 500

   In this test, *10* SNPs are in the segment list. The workspace contains 
   a single annotation. Annotations are all of size *10*, but the overlap of
   SNPs with annotations varies from 0 to 10.



