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

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingGat.png
   :width: 500

   Single continuous workspace.

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

Annotator
---------

.. note::
   If annotator is not installed, the annotator plots will be missing.

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingGat.png
   :width: 500

   Single continuous workspace.

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

Statistics
==========

For 1-sized fragments (i.e. SNPs), the statistics can be checked against
a hypergeometric distribution (sampling without replacement).

Test with a single SNP. Here, there are no issues with multiple hits.

.. figure:: ../test/test_testSingleSNP__main__.TestStatsSNPSampling.png
   :width: 500

In the following test, multiple SNPs are in the segment list, all overlap the
annotations. Hence all results are highly signficant. The size of the annotations 
increases.

.. figure:: ../test/test_testMultipleSNPsFullOverlap__main__.TestStatsSNPSampling.png
   :width: 500

In the following test, multiple SNPs are in the segment list, but with decreasing overlap with the
annotation. Annotations have similar sizes, hence the expected overlap is the same in all tests
and you will expect two point clouds in the top plots. Tests with little observed overlap between 
segments and annotations should not be significant.

.. figure:: ../test/test_testMultipleSNPsPartialOverlap__main__.TestStatsSNPSampling.png
   :width: 500





