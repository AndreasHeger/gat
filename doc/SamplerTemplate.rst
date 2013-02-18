Sampler - %%%
-------------------------------------------

The aim of the plots in this section is to detect
any systematic biases within samplers. The test case
employ different configurations of workspaces and segments
and plot coverage data, segment start points, number of
segments in order to display any biases. The plots 
have four sections. From the top:

1. First plot: Coverage of nucleotides within the workspace by
sampled segments over *n* samples. The default expectation is 
that nucleotides should be equally covered. Gaps within the
workspace have been removed for display purposes. The lines 
show a smoothed coverage curve, the expectation of uniform
coverage with 10% boundaries.

2. Second plot: Segment size distribution. The plot shows summary
statistics for the distribution of segment sizes per sample. The x-axis is 
the sample number. The samples have been sorted by the standard deviation.

3. Third plot: Distribution of segment start and end points. The
x-axis is the same workspace as in the top figure.

4. Fourth plot: The number of segments simulated within each sample.
The plot shows the number of segments input as a dashed line. The
x-axis is the sample number. The samples have been sorted by the
number of segments returned.

Continuous workspaces
+++++++++++++++++++++

.. figure:: ../test/testSingleWorkspace.TestSegmentSampling%%%.png
   :width: 500

   Single continuous workspace.

.. figure:: ../test/testSingleWorkspaceWithOffset.TestSegmentSampling%%%.png
   :width: 500

   Single continuous workspace with an offset. Compare this to the 
   previous plot in order to detect any effects due to workspace
   segments starting at 0.

.. figure:: ../test/testFullWorkspace.TestSegmentSampling%%%.png
   :width: 500

   A single continuous workspace of size 100 containing a single
   segment of size 200. Note the returned segment sizes - as annotator
   will sample until all 100 bases of the workspace are reached, the
   returned segment length can be up to 500 (200 + 100 + 200 ).

.. figure:: ../test/testSmallWorkspace.TestSegmentSampling%%%.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 50.

Segmented workspaces
++++++++++++++++++++

.. figure:: ../test/testSegmentedWorkspaceSmallGap.TestSegmentSampling%%%.png
   :width: 500

   Workspace segmented into 10 segments of size 999 with a single nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspaceLargeGap.TestSegmentSampling%%%.png
   :width: 500

   Workspace segmented into 10 segments of size 900 with a 100 nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspace2x.TestSegmentSampling%%%.png
   :width: 500

   Workspace segmented into 10 segments of size 200 with a 800 nucleotide
   gap between workspaces. In this case, workspace segments are only twice 
   the size of segments.

.. figure:: ../test/testSegmentedWorkspaceSmallGapUnequalSides.TestSegmentSampling%%%.png
   :width: 500

   A segmented workspace of size 100 split at position 50 with a gap of 25. There is 
   a single segment of size 50.

.. figure:: ../test/testSegmentedWorkspaceSmallGapEqualSides.TestSegmentSampling%%%.png
   :width: 500

   A segmented workspace of size 125 split at position 50 with a gap of 5. There is 
   a single segment of size 50.
