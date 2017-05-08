"""test gat implementation."""

import unittest
import os
import numpy

from gat.Engine import AnnotatorResult, IntervalCollection, \
    Samples, SamplesCached, computeFDR, \
    SamplerAnnotator

from gat.SegmentList import SegmentList


class GatTest(unittest.TestCase):

    def shortDescription(self):
        return None


class TestIntervalCollection(GatTest):

    def setUp(self):
        self.a = IntervalCollection("a")
        self.a.add("track1", "contig1",
                   SegmentList(iter=((x, x + 10) for x in range(0, 1000, 100)), normalize=True))
        self.a.add("track1", "contig2",
                   SegmentList(iter=((x, x + 10) for x in range(0, 1000000, 100)), normalize=True))
        self.a.add("track2", "contig1",
                   SegmentList(iter=((x, x + 10) for x in range(1000, 2000, 100)), normalize=True))

    def testSaveLoad(self):

        fn = "tmp_testSaveLoad.bed"
        outfile = open(fn, "w")
        self.a.save(outfile)
        outfile.close()

        b = IntervalCollection("b")
        b.load(fn)

        self.assertEqual(len(self.a), len(b))
        self.assertEqual(self.a.tracks, b.tracks)
        self.assertEqual(self.a.sum(), b.sum())
        for t in b.tracks:
            self.assertEqual(sorted(self.a[t].keys()), sorted(b[t].keys()))
            for x in list(self.a[t].keys()):
                self.assertEqual(self.a[t][x], b[t][x])

    def testSharing(self):

        aa = self.a.clone()
        aa.share()


class TestToFromIsochores(GatTest):

    def setUp(self):

        a = IntervalCollection("a")

        # every ten at
        a.add("track1", "contig1",
              SegmentList(iter=((x, x + 250) for x in range(0, 2000, 500)), normalize=True))

        a.add("track1", "contig2",
              SegmentList(iter=((x, x + 250) for x in range(0, 2000, 500)), normalize=True))

        self.a = a

    def check(self, a, orig):

        self.assertEqual(a.tracks, orig.tracks)
        self.assertEqual(sorted(a["track1"].keys()),
                         sorted(['contig2', 'contig1']))

        for contig in ("contig1", "contig2"):
            self.assertEqual(list(a["track1"][contig]),
                             list(orig["track1"][contig]))

    def testNoIsochores(self):

        a = self.a.clone()
        orig = a.clone()
        a.fromIsochores()
        self.check(a, orig)

    def testToFromIsochores(self):

        a = self.a.clone()
        orig = a.clone()

        isochores = IntervalCollection("isochores")

        # covering isochores of size 100
        isochores.add("highGC", "contig1",
                      SegmentList(iter=((x, x + 100) for x in range(0, 10000, 200)), normalize=True))

        isochores.add("lowGC", "contig1",
                      SegmentList(iter=((x, x + 100) for x in range(100, 10000, 200)), normalize=True))

        isochores.add("highGC", "contig2",
                      SegmentList(iter=((x, x + 100) for x in range(0, 10000, 200)), normalize=True))

        isochores.add("lowGC", "contig2",
                      SegmentList(iter=((x, x + 100) for x in range(100, 10000, 200)), normalize=True))

        a.toIsochores(isochores)
        self.assertEqual(a.tracks, orig.tracks)
        self.assertEqual(sorted(a["track1"].keys()),
                         sorted(['contig2.highGC', 'contig1.highGC', 'contig2.lowGC', 'contig1.lowGC']))

        a.fromIsochores()

        self.check(a, orig)


class TestPValue(GatTest):
    '''test if pvalue computation is correct.'''

    def testPValue1(self):

        observed = 0.332640195285
        values = [0.3593727449353678, 0.24446041723385858, 0.11321078358680142, 0.28500665546177717, 0.017634423032620888, 0.47144573882791929, 0.20295266762886535, 0.24374906675401431, 0.12536987767373536, 0.36647407597049514, 0.1317950839045143, 0.32036858313479905, 0.2131875486832529, 0.18211958887292382, 0.4382662865088186, 0.12487068923091568, 0.38895983423268921, 0.43156050120631062, 0.18784825518278428, 0.23958644581530344, 0.16386055449534453, 0.42697777787951602, 0.07748674945294963, 0.47881248869131277, 0.37267534771232319, 0.8083924735050152, 0.29179189019925428, 0.29802029242077777, 0.2054027587360118, 0.10766996738143179, 0.39134998593956405, 0.36412616130029274, 0.37015995608450686, 0.61246049427537563, 0.59897086243095388, 0.20718454055122912, 0.14334918487088333, 0.42189815231899974, 0.21738749430714899, 0.39304902005163428, 0.50261637732761, 0.20759334134444557, 0.21005124432686503, 0.31027042275886835, 0.71335371670327341, 1.4192781030245714, 0.50672517580861098, 0.18067653694488042, 0.85952730574991043, 0.19249388587333111,
                  0.18826477050167958, 0.22742885130411533, 0.24125995809534906, 0.045750800392306591, 0.78242626285998884, 0.20614461737324383, 0.56783904985512668, 0.33500622312674566, 0.043533317315170454, 0.27874382197104552, 0.3685525858770754, 0.1751812314517863, 0.2532293526642409, 0.15785104775566922, 0.2390711833181299, 0.42911409505776471, 0.16819203200742916, 0.40372196988518594, 0.43241512178368696, 0.30424021778439686, 0.19085162018033855, 0.58462246847853661, 0.631050399423982, 0.30137309454374051, 0.27565096287611918, 0.33033553618821287, 0.47665164689105288, 0.34084703029218633, 0.27978844627773986, 0.010536582324145049, 0.050935298127348511, 0.23536721808983668, 0.22364077067346355, 0.31704429093519465, 1.0296141286403104, 0.38123158028252929, 0.27538594123938104, 0.81446088474774558, 0.2660327021825486, 0.195234318277725, 0.462999083371401, 1.0587870384537303, 0.40260375543813692, 0.39471961997665139, 0.29845505700189406, 1.0259474557694457, 0.52111381852233729, 0.29182304834212835, 0.34045457181768657, 0.20417518807825608]

        result = AnnotatorResult(
            "test", "test", "na", observed, values, reference=None, pseudo_count=0)

        self.assertEqual(result.pvalue, 0.57)


class TestSamples(GatTest):

    nsamples = 1000
    ntracks = 50
    nsegments = 10000
    nisochores = 20

    def testDelete(self):
        '''test to track down a memory leak.'''
        # from guppy import hpy
        # hp = hpy()
        # hp.setrelheap()
        return

        samples = Samples()

        for track in range(self.ntracks):
            track_id = str(track)
            # print track_id
            for sample in range(self.nsamples):
                sample_id = str(sample)
                for isochore in range(self.nisochores):
                    isochore_id = str(isochore)
                    r = SegmentList(allocate=self.nsegments)
                    samples.add(track_id, sample_id, isochore_id, r)

            del samples[track_id]
            # print len(samples)
            # h = hp.heap()

            # print h


class TestCaching(GatTest):

    sample_size = 10
    workspace_size = 1000

    def testCaching(self):
        return
        workspaces, segments, annotations = \
            IntervalCollection( "workspace" ), \
            IntervalCollection( "segment" ), \
            IntervalCollection("annotation")

        workspaces.add("default", "chr1", SegmentList(iter=[(0, self.workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        segments.add("default", "chr1", SegmentList(iter=[(0, 1), ],
                                                    normalize=True))

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(1, 100, 2):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(0, y), ],
                                        normalize=True))

        workspace_size = workspace["chr1"].sum()

        sampler = SamplerAnnotator(bucket_size=1, nbuckets=self.workspace_size)

        if os.path.exists("test.cache"):
            os.remove("test.cache")

        outsamples = SamplesCached("test.cache")
        saved_samples = {}

        for track in segments.tracks:
            segs = segments[track]
            for x in range(self.sample_size):
                for isochore in list(segs.keys()):
                    r = sampler.sample(segs[isochore], workspace[isochore])
                    saved_samples[(track, x, isochore)] = r
                    outsamples.add(track, x, isochore, r)

        del outsamples

        insamples = SamplesCached("test.cache")

        for track in segments.tracks:
            segs = segments[track]
            for x in range(self.sample_size):
                for isochore in list(segs.keys()):
                    insamples.load(track, x, isochore)

        # compare
        for track in segments.tracks:
            segs = segments[track]
            for x in range(self.sample_size):
                for isochore in list(segs.keys()):
                    self.assertEqual(saved_samples[(track, x, isochore)].asList(),
                                     insamples[track][x][isochore].asList())

        if os.path.exists("test.cache"):
            os.remove("test.cache")
        if os.path.exists("test.cache.idx"):
            os.remove("test.cache.idx")


class TestStats(GatTest):

    ntracks = 10  # 17
    nannotations = 10  # 90
    nsamples = 1000

    def testPValueComputation(self):

        annotation_size = 100
        workspace_size = 1000
        segment_size = 10

        l = 10

        for y in range(1, l):

            samples = [1] * y + [0] * (l - y)

            for x, s in enumerate(samples):

                g = AnnotatorResult("track",
                                    "samples",
                                    'counter',
                                    s,
                                    samples)
                self.assertEqual(
                    g.isSampleSignificantAtPvalue(x, g.pvalue), True)

                t = 0
                for z, s2 in enumerate(samples):
                    t += g.isSampleSignificantAtPvalue(z, g.pvalue)
                fpr = float(t) / l

                # == should work, but there is a problem
                # for pvalue = 0.5
                self.assertTrue(fpr >= g.pvalue - 0.0001,
                             "fpr (%f) != pvalue (%f): y=%i, x=%i, s=%i, t=%i" %
                             (fpr, g.pvalue, y, x, s, t))

    def testPValueComputation2(self):

        samples = [0] * 66 + [1] * 2 + [2] * 20 + [3] * \
            1 + [4] * 6 + [6] * 2 + [8] * 2 + [16] * 1

        obs = 16

        g = AnnotatorResult("track", "samples",
                            'counter',
                            obs,
                            samples)

        self.assertEqual(g.pvalue, 0.01)

    def testStats(self):

        annotation_size = 100
        workspace_size = 1000
        segment_size = 10

        observed = numpy.random.hypergeometric(annotation_size,
                                               workspace_size - annotation_size,
                                               segment_size,
                                               self.ntracks * self.nannotations)

        results = []
        x = 0
        for track in range(self.ntracks):
            for annotation in range(self.nannotations):
                samples = numpy.random.hypergeometric(annotation_size,
                                                      workspace_size - annotation_size,
                                                      segment_size,
                                                      self.nsamples)

                samples.sort()
                # print "obs", observed[x], "sample=", samples

                results.append(AnnotatorResult(str(track),
                                               str(annotation),
                                               'counter',
                                               observed[x],
                                               samples))

                x += 1

        computeFDR(results)
        for r in results:
            self.assertTrue(r.qvalue > 0.5, "%f" % r.qvalue)

if __name__ == '__main__':
    unittest.main()
