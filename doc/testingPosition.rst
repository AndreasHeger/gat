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

