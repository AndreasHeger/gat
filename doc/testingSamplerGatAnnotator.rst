SamplerAnnotator
----------------

Continuous workspaces
+++++++++++++++++++++

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingGat.png
   :width: 500

   Single continuous workspace.

.. figure:: ../test/testFullWorkspace.TestSegmentSamplingGat.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 200. Note the returned segment sizes - as annotator
   will sample until all 100 bases of the workspace are reached, the
   returned segment length can be up to 500 (200 + 100 + 200 ).

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
