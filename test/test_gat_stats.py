'''test gat functionality
'''

import unittest
import gat


class TestSNPs(unittest.TestCase):

    def check(self):
        pass

    def testIntervalsPartialOverlap(self):
        '''test with intervals with 
        increasing amount of overlap.

        '''
        return
        workspaces, segments, annotations = \
            gat.IntervalCollection("workspace"), \
            gat.IntervalCollection("segment"), \
            gat.IntervalCollection("annotation")

        workspace_size = 1000

        size = 100

        # workspace of size 1000000
        workspaces.add("default", "chr1",
                       gat.SegmentList(iter=[(0, workspace_size), ],
                                       normalize=True))
        workspace = workspaces["default"]

        # segment of size 10
        segments.add("default", "chr1",
                     gat.SegmentList(iter=[(0, size), ],
                                     normalize=True))

        # annotations: a collection of segments.
        # overlap increases
        annotations.add("full", "chr1",
                        gat.SegmentList(iter=[(y, size + y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)


if __name__ == '__main__':
    unittest.main()
