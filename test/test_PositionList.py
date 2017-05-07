"""test gat implementation."""

import unittest
import random
import os
import numpy
import pickle

from gat.PositionList import PositionList
from gat.SegmentList import SegmentList


class GatTest(unittest.TestCase):

    def shortDescription(self):
        return None


class TestPositionList(GatTest):
    '''test segment list implementation.'''

    def testCreateAndClear(self):
        s = PositionList()
        self.assertEqual(0, len(s))
        s.add(100)
        self.assertEqual(1, len(s))
        s.clear()
        self.assertEqual(0, len(s))

    def testNormalize(self):
        '''non-overlapping segments.'''

        ss = list(range(0, 1000, 100))
        random.shuffle(ss)
        s = PositionList()
        for pos in ss:
            s.add(pos)
        self.assertEqual(len(s), 10)
        s.normalize()
        self.assertEqual(len(s), 10)
        self.assertEqual(s.sum(), 10)
        for pos in ss:
            s.add(pos)
        self.assertEqual(len(s), 20)
        s.normalize()
        self.assertEqual(len(s), 10)

    def testFromSegments(self):
        l = [(489, 589),
             (1966, 2066),
             (2786, 2886),
             (3889, 3972),
             (3998, 4098),
             (6441, 6541),
             (6937, 7054),
             (7392, 7492),
             (8154, 8254),
             (9046, 9146)]
        ss = SegmentList(iter=l)
        pp = PositionList()
        pp.fromSegmentList(ss)
        self.assertEqual(len(ss), len(pp))

    def testOverlapWithRange(self):

        # single point per position
        pp = PositionList(iter=list(range(0, 1000, 100)), sort=True)

        for x in range(0, 1050, 10):
            if x % 100 == 0:
                if x == 1000:
                    self.assertEqual(0, pp.overlapWithRange(x, x+10))
                else:
                    self.assertEqual(1, pp.overlapWithRange(x, x+10))
            else:
                self.assertEqual(0, pp.overlapWithRange(x, x+10))

        # two points per positions
        pp = PositionList(iter=list(range(0, 1000, 100)) + list(range(0, 1000, 100)),
                          sort=True)

        for x in range(0, 1050, 10):
            if x % 100 == 0:
                if x == 1000:
                    self.assertEqual(0, pp.overlapWithRange(x, x+10))
                else:
                    self.assertEqual(2, pp.overlapWithRange(x, x+10))
            else:
                self.assertEqual(0, pp.overlapWithRange(x, x+10))

    def testOverlapWithSegments(self):

        # single point per position
        pp = PositionList(iter=list(range(0, 1000, 100)), sort=True)

        for o in range(0, 200, 10):
            ss = SegmentList(
                iter=[(x, x + 1) for x in range(0 + o, 1000 + o, 100)],
                normalize=True)
            if o % 100 == 0:
                if o == 100:
                    self.assertEqual(9,
                                     pp.intersectionWithSegments(ss))
                else:
                    self.assertEqual(10,
                                     pp.intersectionWithSegments(ss))
            else:
                self.assertEqual(0,
                                 pp.intersectionWithSegments(ss))

    def testIntersect(self):

        for o in range(0, 200, 10):
            # single point per position
            pp = PositionList(iter=list(range(0, 1000, 100)), sort=True)

            ss = SegmentList(
                iter=[(x, x + 1) for x in range(0 + o, 1000 + o, 100)],
                normalize=True)
            pp.intersect(ss)
            if o % 100 == 0:
                if o == 100:
                    self.assertEqual(9, len(pp))
                else:
                    self.assertEqual(10, len(pp))
            else:
                self.assertEqual(0, len(pp))

if __name__ == '__main__':
    unittest.main()
