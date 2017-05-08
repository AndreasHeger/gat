"""test gat implementation."""

import unittest
import random
import pickle

from gat.SegmentList import SegmentList


class GatTest(unittest.TestCase):

    def shortDescription(self):
        return None


class TestSegmentList(GatTest):
    '''test segment list implementation.'''

    def testCreateAndClear(self):
        s = SegmentList()
        self.assertEqual(0, len(s))
        s.add(0, 100)
        self.assertEqual(1, len(s))
        s.clear()
        self.assertEqual(0, len(s))

    def testNormalize1(self):
        '''non-overlapping segments.'''

        ss = [(x, x + 10) for x in range(0, 1000, 100)]
        random.shuffle(ss)
        s = SegmentList()
        for start, end in ss:
            s.add(start, end)
        s.normalize()

        self.assertEqual(len(s), 10)
        self.assertEqual(s.sum(), 100)
        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalize1b(self):
        '''non-overlapping segments.'''

        ss = [(x, x + 10) for x in range(100, 1100, 100)]
        random.shuffle(ss)
        s = SegmentList()
        for start, end in ss:
            s.add(start, end)
        s.normalize()
        self.assertEqual(len(s), 10)
        self.assertEqual(s.sum(), 100)
        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalizeEmpty(self):
        '''non-overlapping segments.'''

        s = SegmentList()
        self.assertEqual(len(s), 0)
        s.normalize()
        self.assertEqual(len(s), 0)
        self.assertEqual(s.isNormalized, 1)
        s2 = SegmentList()
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalizeEmptySegment(self):
        s = SegmentList(iter=[(0, 0), ])
        s.normalize()
        self.assertEqual(s.isNormalized, 1)
        self.assertEqual(len(s), 0)

        s = SegmentList(iter=[(0, 0), (0, 0)])
        s.normalize()
        self.assertEqual(s.isNormalized, 1)
        self.assertEqual(len(s), 0)

        ss = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
              (0, 5), (0, 6), (0, 7), (0, 8), (0, 9)]
        s = SegmentList(iter=ss)
        s.normalize()
        self.assertEqual(s.isNormalized, 1)
        self.assertEqual(len(s), 1)

        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalize2(self):
        '''overlapping segments.'''

        ss = [(x, x + 1000) for x in range(0, 1000, 100)]
        random.shuffle(ss)
        s = SegmentList()
        for start, end in ss:
            s.add(start, end)
        s.normalize()
        self.assertEqual(len(s), 1)
        self.assertEqual(s.sum(), 1900)

    def testNormalize3(self):
        '''non-overlapping but adjacent segments.'''

        ss = [(x, x + 100) for x in range(0, 1000, 100)]
        random.shuffle(ss)
        s = SegmentList()
        for start, end in ss:
            s.add(start, end)
        s.normalize()

        self.assertEqual(len(s), 10)
        self.assertEqual(s.sum(), 1000)
        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalize4(self):
        # test multiple interleaved segments
        ss = [(x, x + 100) for x in range(0, 1000, 10)]
        s = SegmentList()
        for start, end in ss:
            s.add(start, end)
        s.normalize()
        self.assertEqual(len(s), 1)
        self.assertEqual(s.sum(), 1090)
        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testNormalize5(self):
        ss = [(489, 589), (1966, 2066), (2786, 2886), (0, 0), (3889, 3972), (3998, 4098),
              (6441, 6541), (6937, 7054), (7392, 7492), (8154, 8254), (9046, 9146)]
        s = SegmentList(iter=ss)
        s.normalize()
        self.assertEqual(len(s), len([x for x in ss if x[1] - x[0] > 0]))
        self.assertEqual(s.sum(), 1000)
        s2 = SegmentList(iter=ss)
        s2.merge(-1)
        self.assertEqual(s, s2)

    def testExtend(self):

        s1 = SegmentList(iter=[(x, x + 100) for x in range(0, 1000, 100)])
        s2 = SegmentList(iter=[(x, x + 100) for x in range(2000, 3000, 100)])
        s1.extend(s2)
        self.assertEqual(s1.sum(), s2.sum() * 2)
        self.assertEqual(len(s1), len(s2) * 2)

    def testTrim(self):
        '''test trimming over full range of insertion points and deletions.'''

        for point in range(0, 1000):
            for size in range(0, 300):
                ss = [(x, x + 100) for x in range(0, 1000, 100)]
                s = SegmentList(iter=ss,
                                normalize=True)
                orig = s.sum()
                s.trim(point, size)
                self.assertEqual(orig - size, s.sum(),
                                 "trimming error at %i:%i: expected %i, got %i, %s" %
                                 (point, size,
                                  orig - size, s.sum(), str(s)))

    def testInsertionPoint(self):
        '''check insertion point for normalized segment lists'''
        ss = [(x, x + 10) for x in range(0, 100, 10)]
        s = SegmentList(iter=ss,
                        normalize=True)

        for point in range(0, 100):
            p = s.getInsertionPoint(point, point + 1)
            self.assertEqual(p, point // 10)

        ss = [(x, x + 10) for x in range(0, 100, 20)]
        s = SegmentList(iter=ss,
                        normalize=True)

        for point in range(0, 100):
            p = s.getInsertionPoint(point, point + 1)
            if point >= 90:
                self.assertEqual(p, len(s))
            else:
                self.assertEqual(p, point // 20)

        ss = [(x, x + 10) for x in range(10, 100, 20)]
        s = SegmentList(iter=ss,
                        normalize=True)

        for point in range(0, 100):
            p = s.getInsertionPoint(point, point + 1)
            self.assertEqual(p, (point - 10) // 20)

    def testInsertionPointNonNormalized(self):
        '''check insertion point for unnormalized segment lists.'''
        ss = [(x, x + 20) for x in range(0, 100, 10)]
        s = SegmentList(iter=ss,
                        normalize=False)

        for point in range(0, 100):
            self.assertRaises(
                AssertionError, s.getInsertionPoint, point, point + 1)

    def testMergeAdjacent(self):
        ss = [(x, x + 100) for x in range(0, 1000, 100)]
        random.shuffle(ss)
        s = SegmentList(iter=ss)
        s.merge(0)
        self.assertEqual(len(s), 1)
        self.assertEqual(s.sum(), 1000)

    def testMergeNeighbours(self):

        for y in range(0, 5):
            ss = [(x, x + 100 - y) for x in range(0, 1000, 100)]
            random.shuffle(ss)
            for x in range(0, y + 1):
                s = SegmentList(iter=ss)
                s.merge(x)
                if x < y:
                    self.assertEqual(len(s), 10)
                    self.assertEqual(s.sum(), 1000 - 10 * y)
                else:
                    self.assertEqual(len(s), 1)
                    self.assertEqual(s.sum(), 1000 - y)

    def testExpand(self):
        ss = [(x, x + 10) for x in range(0, 1000, 100)]
        s = SegmentList(iter=ss, normalize=True)
        s.expand_segments(2.0)
        self.assertEqual(s.sum(), 195)
        self.assertEqual(len(s), 10)

    def testLargest(self):
        ss = [(x, x + x / 10) for x in range(0, 1000, 100)]
        s = SegmentList(iter=ss, normalize=True)
        p = s.largest()
        self.assertEqual(p['end'] - p['start'], 90)

    def testGetFilledSegmentsFromStart(self):
        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)
        for x in range(0, 120, 5):
            f = s.getFilledSegmentsFromStart(x, 20)
            self.assertEqual(f.sum(), 20)
            if x >= 110:
                self.assertEqual(f.min(), s.min())
                self.assertEqual(f.max(), 30)
            elif x > 80:
                self.assertEqual(f.min(), s.min())
                self.assertEqual(f.max(), s.max())
            else:
                if x in (0, 20, 40, 60, 80):
                    self.assertEqual(f.min(), x)
                    self.assertEqual(f.max(), x + 30)
                elif x in (10, 30, 50, 70):
                    self.assertEqual(f.min(), x + 10)
                    self.assertEqual(f.max(), x + 40)
                elif x in (5, 25, 45, 65):
                    self.assertEqual(f.min(), x)
                    self.assertEqual(f.max(), x + 40)
                elif x in (15, 35, 55, 75):
                    self.assertEqual(f.min(), x + 5)
                    self.assertEqual(f.max(), x + 35)
            self.assertEqual(
                s.sum(), s.getFilledSegmentsFromStart(x, 100).sum())

    def testGetFilledSegmentsFromEnd(self):
        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)
        for x in range(0, 120, 5):
            f = s.getFilledSegmentsFromEnd(x, 20)
            self.assertEqual(f.sum(), 20)
            if x == 0:
                self.assertEqual(f.max(), s.max())
                self.assertEqual(f.min(), 80)
            elif x < 30:
                self.assertEqual(f.max(), s.max())
                self.assertEqual(f.min(), s.min())
            else:
                if x in (40, 60, 80, 100):
                    self.assertEqual(f.min(), x - 40)
                    self.assertEqual(f.max(), x - 10)
                elif x in (30, 50, 70, 90):
                    self.assertEqual(f.min(), x - 30)
                    self.assertEqual(f.max(), x)
                elif x in (45, 65, 85):
                    self.assertEqual(f.min(), x - 40)
                    self.assertEqual(f.max(), x)
                elif x in (35, 55, 75, 95):
                    self.assertEqual(f.min(), x - 35)
                    self.assertEqual(f.max(), x - 5)

            self.assertEqual(s.sum(), s.getFilledSegmentsFromEnd(x, 100).sum())

    def testPickling(self):
        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)

        b = pickle.loads(pickle.dumps(s))

        self.assertEqual(s, b)

    def testSharing(self):

        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)
        s.share("/testshare")
        n = SegmentList(share=s)
        self.assertEqual(s, n)

    def testPickledSharing(self):

        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)
        s.share("/testshare")
        b = pickle.loads(pickle.dumps(s))
        self.assertEqual(s, b)

    def testUnsharing(self):

        ss = [(x, x + 10) for x in range(0, 120, 20)]
        s = SegmentList(iter=ss, normalize=True)
        s.share("/testshare")
        s.unshare()
        s.share("/testshare")
        s.unshare()


class TestSegmentListOverlap(GatTest):

    def setUp(self):
        self.a = SegmentList(iter=((x, x + 10)
                                   for x in range(0, 1000, 100)), normalize = True)

    def testOverlapFull(self):
        self.assertEqual(self.a.overlapWithRange(0, 1000), self.a.sum())

    def testOverlapHalf(self):
        self.assertEqual(self.a.overlapWithRange(0, 500), self.a.sum() / 2)
        self.assertEqual(self.a.overlapWithRange(500, 1000), self.a.sum() / 2)

    def testOverlapAfter(self):
        self.assertEqual(self.a.overlapWithRange(900, 910), 10)
        self.assertEqual(self.a.overlapWithRange(905, 915), 5)

    def testOverlapNone(self):
        self.assertEqual(self.a.overlapWithRange(1000, 2000), 0)
        self.assertEqual(self.a.overlapWithRange(2000, 3000), 0)

    def testOverlapOutOfRange(self):
        self.assertRaises(OverflowError, self.a.overlapWithRange, -100, 5)
        self.assertRaises(OverflowError, self.a.overlapWithRange, -100, -50)
        self.assertEqual(self.a.overlapWithRange(905, 1100), 5)
        self.assertEqual(self.a.overlapWithRange(1000, 1100), 0)

    def testOverlapAll(self):
        for x in range(0, 1000, 100):
            for y in range(0, 10):
                self.assertEqual(self.a.overlapWithRange(x + y, x + y + 1), 1,
                                 "no overlap failure at %i: %i" % (x + y, self.a.overlapWithRange(x + y, x + y + 1)))
            for y in range(10, 100):
                self.assertEqual(self.a.overlapWithRange(
                    x + y, x + y + 1), 0, "overlap failure at %i" % (x + y))


class TestSegmentListIntersection(GatTest):

    def setUp(self):
        self.a = SegmentList(
            iter=((x, x + 10) for x in range(0, 1000, 100)), normalize=True)

    def testIntersectionFull(self):
        b = SegmentList(iter=[(0, 1000)], normalize=True)
        b.intersect(self.a)
        self.assertEqual(b.asList(), self.a.asList())

    def testIntersectionSelf(self):
        self.a.intersect(self.a)
        self.assertEqual(self.a.asList(), self.a.asList())

    def testIntersectionCopy(self):
        b = SegmentList(clone=self.a)
        b.intersect(self.a)
        self.assertEqual(b.asList(), self.a.asList())

    def testNoIntersection(self):
        b = SegmentList(iter=((x, x + 10)
                              for x in range(10, 1000, 100)), normalize = True)
        b.intersect(self.a)
        self.assertEqual(b.asList(), [])
        self.assertEqual(b.isEmpty, True)

    def testPartialIntersection(self):
        b = SegmentList(iter=((x, x + 10)
                              for x in range(5, 1000, 100)), normalize = True)
        b.intersect(self.a)
        self.assertEqual(len(b), len(self.a))
        self.assertEqual(b.sum(), self.a.sum() / 2)

    def testOverlap(self):
        '''test if number of segments intersection is correct.'''

        b = SegmentList(iter=((x, x + 10)
                              for x in range(5, 1000, 100)), normalize = True)
        self.assertEqual(self.a.intersectionWithSegments(b), len(b))
        self.assertEqual(b.intersectionWithSegments(self.a), len(b))

        # no intersection
        b = SegmentList(iter=((x, x + 10)
                              for x in range(10, 1000, 100)), normalize = True)
        self.assertEqual(self.a.intersectionWithSegments(b), 0)
        self.assertEqual(b.intersectionWithSegments(self.a), 0)

        # double the number of segments in b
        b = SegmentList(iter=[(x, x + 5) for x in range(0, 1000, 100)] +
                        [(x + 5, x + 10) for x in range(0, 1000, 100)],
                        normalize=True)
        self.assertEqual(self.a.intersectionWithSegments(b), 10)
        self.assertEqual(b.intersectionWithSegments(self.a), 20)

    def testFilter(self):

        b = SegmentList(iter=((x, x + 5)
                              for x in range(500, 2000, 100)), normalize=True)
        b.filter(self.a)
        self.assertEqual(b.asList(), [
            (500, 505), (600, 605), (700, 705), (800, 805),
            (900, 905)])

        b = SegmentList(iter=((0, 56), ))
        c = SegmentList(iter=[(0, 50), (75, 125)])
        b.filter(c)
        self.assertEqual(b.asList(), [(0, 56)])

        b = SegmentList(iter=((0, 56), ))
        c = SegmentList(iter=[(0, 10)])
        b.filter(c)
        self.assertEqual(b.asList(), [(0, 56)])


class TestSegmentListSubtract(GatTest):

    def setUp(self):
        self.a = SegmentList(
            iter=((x, x + 10) for x in range(0, 1000, 100)), normalize=True)

    def testCompleteOverlap(self):
        b = SegmentList(iter=[(0, 1000)], normalize=True)
        b.subtract(self.a)
        c = [(10, 100), (110, 200), (210, 300), (310, 400), (410, 500),
             (510, 600), (610, 700), (710, 800), (810, 900), (910, 1000)]
        self.assertEqual(b.asList(), c)

    def testFullSubtraction(self):
        b = SegmentList(iter=[(0, 1000)], normalize=True)
        self.a.subtract(b)
        self.assertEqual(len(self.a), 0)

    def testSelfSubtraction(self):
        self.a.subtract(self.a)
        self.assertEqual(len(self.a), 0)

    def testSameSubtraction(self):
        b = SegmentList(clone=self.a)
        b.subtract(self.a)
        self.assertEqual(len(b), 0)

    def testOverlap(self):

        b = SegmentList(iter=((x, x + 10)
                              for x in range(5, 1000, 100)), normalize = True)
        b.subtract(self.a)
        c = [(10, 15), (110, 115), (210, 215), (310, 315), (410, 415),
             (510, 515), (610, 615), (710, 715), (810, 815), (910, 915)]
        self.assertEqual(b.asList(), c)

    def testSingleSegmentSubtraction(self):
        a = SegmentList(iter=[(0, 12000)], normalize=True)
        b = SegmentList(iter=[(0, 10000)], normalize=True)
        a.subtract(b)
        self.assertEqual(a.asList(), [(10000, 12000)])


if __name__ == '__main__':
    unittest.main()
