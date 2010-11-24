============
Benchmarking
============

This section contains quality control plots from the unit testing.

Position sampling
=================

The :class:`TestSamplerPosition` tests if the position sampling works
as expected. In particular, look out for edge effects.

.. figure:: ../test/test_testMultipleWorkspaces__main__.TestSamplerPosition.png
   :width: 500

   Multiple work spaces. 10 workspaces of size 100, spaced every 1000 nucleotides

.. figure:: ../test/test_testPositionSamplingSingleWorkspace__main__.TestSamplerPosition.png
   :width: 500

   A single work space.

.. figure:: ../test/test_testPositionSamplingSplitWorkspace__main__.TestSamplerPosition.png
   :width: 500

   A workspace split in the middle without a gap.

.. figure:: ../test/test_testPositionSamplingSplitWorkspace2__main__.TestSamplerPosition.png
   :width: 500

   A workspace split in the middle with a gap in between.

.. figure:: ../test/test_testSNPPositionSampling__main__.TestSamplerPosition.png
   :width: 500

   10 workspaces of size 100, segment size of 1 (SNP).

Sampling
========

The following plots benchmark the segment sampling behaviour of GAT and Annotator.
The tests below sample from a segment list *10* segments with each segment of size *100*.


GAT
---

Continuous workspaces
+++++++++++++++++++++

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingGat.png
   :width: 500

   Single continuous workspace.

.. figure:: ../test/testFullWorkspace.TestSegmentSamplingGat.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 200.

.. figure:: ../test/testSmallWorkspace.TestSegmentSamplingGat.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 50.

Segmented workspaces
++++++++++++++++++++

.. figure:: ../test/testSegmentedWorkspaceSmallGap.TestSegmentSamplingGat.png
   :width: 500

   Workspace segmented into 10 segments of size 999 with a single nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspaceLargeGap.TestSegmentSamplingGat.png
   :width: 500

   Workspace segmented into 10 segments of size 900 with a 100 nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspace2x.TestSegmentSamplingGat.png
   :width: 500

   Workspace segmented into 10 segments of size 200 with a 800 nucleotide
   gap between workspaces. In this case, workspace segments are only twice 
   the size of segments.

.. figure:: ../test/testSegmentedWorkspaceSmallGapUnequalSides.TestSegmentSamplingGat.png
   :width: 500

   A segmented workspace of size 100 split at position 50 with a gap of 25. There is 
   a single segment of size 50.

.. figure:: ../test/testSegmentedWorkspaceSmallGapEqualSides.TestSegmentSamplingGat.png
   :width: 500

   A segmented workspace of size 125 split at position 50 with a gap of 5. There is 
   a single segment of size 50.

Annotator
---------

.. note::
   If annotator is not installed, the annotator plots will be missing.

Continuous workspaces
+++++++++++++++++++++

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingTheAnnotator.png
   :width: 500

   Single continuous workspace.

.. figure:: ../test/testFullWorkspace.TestSegmentSamplingTheAnnotator.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 200.

.. figure:: ../test/testSmallWorkspace.TestSegmentSamplingTheAnnotator.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 50.

Segmented workspaces
++++++++++++++++++++

.. figure:: ../test/testSegmentedWorkspaceSmallGap.TestSegmentSamplingTheAnnotator.png
   :width: 500

   Workspace segmented into 10 segments of size 999 with a single nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspaceLargeGap.TestSegmentSamplingTheAnnotator.png
   :width: 500

   Workspace segmented into 10 segments of size 900 with a 100 nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspace2x.TestSegmentSamplingTheAnnotator.png
   :width: 500

   Workspace segmented into 10 segments of size 200 with a 800 nucleotide
   gap between workspaces. In this case, workspace segments are only twice 
   the size of segments.

.. figure:: ../test/testSegmentedWorkspaceSmallGapUnequalSides.TestSegmentSamplingTheAnnotator.png
   :width: 500

   A segmented workspace of size 100 split at position 50 with a gap of 25. There is 
   a single segment of size 50.

.. figure:: ../test/testSegmentedWorkspaceSmallGapEqualSides.TestSegmentSamplingTheAnnotator.png
   :width: 500

   A segmented workspace of size 125 split at position 50 with a gap of 5. There is 
   a single segment of size 50.

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



