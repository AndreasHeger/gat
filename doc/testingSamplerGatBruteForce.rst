SamplerBruteForce
-----------------

Continuous workspaces
+++++++++++++++++++++

.. figure:: ../test/testSingleWorkspace.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   Single continuous workspace.

.. figure:: ../test/testFullWorkspace.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 200.

.. figure:: ../test/testSmallWorkspace.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   A single continuous workspace of size 100. Samples contain a single
   segment of size 50.

.. figure:: ../test/testTinyWorkspace.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   A single continuous workspace of size 12. Samples contain a single
   segment of size 4.

Note the bias within the tiny workspace::

       0123456789012 workspace
       |           |
    ####        ####
     ####        ####
      ####        ####
       ####
	####
	 ####
	  ####
	   ####
	    ####
	     ####
	      ####
	       ####
    ####        ####
     ####        ####
      ####        ####
       |           |
       0123456789012

       766544444567 counts

Sampling starting position of segments uniformly gives rise to a density bias.
The bias is large here as the chance of hitting a workspace boundary is high
due to the small workspace size with respect to the segment size.

Segmented workspaces
++++++++++++++++++++

.. figure:: ../test/testSegmentedWorkspaceSmallGap.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   Workspace segmented into 10 segments of size 999 with a single nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspaceLargeGap.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   Workspace segmented into 10 segments of size 900 with a 100 nucleotide
   gap between workspaces.

.. figure:: ../test/testSegmentedWorkspace2x.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   Workspace segmented into 10 segments of size 200 with a 800 nucleotide
   gap between workspaces. In this case, workspace segments are only twice 
   the size of segments.

.. figure:: ../test/testSegmentedWorkspaceSmallGapUnequalSides.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   A segmented workspace of size 100 split at position 50 with a gap of 25. There is 
   a single segment of size 50.

.. figure:: ../test/testSegmentedWorkspaceSmallGapEqualSides.TestSegmentSamplingSamplerBruteForce.png
   :width: 500

   A segmented workspace of size 125 split at position 50 with a gap of 5. There is 
   a single segment of size 50.

