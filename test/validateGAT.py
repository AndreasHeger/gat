'''GAT needs to

* account for length biases
* account of sequence composition (G+C content)

Default case
============

Annotations: There are 10 sets of 10 genes of incremental size.

Segments: Segments are uniformly distributed over the complete workspace.

Workspace: Complete

Expected result: No annotations will be reported as enriched.

Ascertainment bias
==================

Not all genomic regions are always amenable for segments. Setting the size of the workspace
appropriately can correct this bias.

testAscertainmentBiasFail1: The workspace is too large
* Annotations: There are 10 sets of 10 annotations of incremental size.
* Segments: Segments are uniformly distributed over complete workspace
* Workspace: Twice the size of the complete workspace
* Expected result: All annotations will be reported as enriched.

testAscertainmentBiasFail2: Segments only in annotated workspace

* Annotations: There are 10 sets of 10 annotations of incremental size.
* Segments: Segments are uniformly distributed over annotated workspace.
* Workspace: Complete
* Expected result: All annotations will be reported as enriched.

testAscertainmentBiasPass1: Workspace restricted to annotations.
* Workspace: Restricted to annotations
* Expected result: No annotations will be reported as enriched.

Length bias
============

Longer annotations have a higher chance of being sampled by segments. When computing 
enrichment, a test based on overlap counts (such as in a gene list analysis) will call longer 
annotations to be enriched.

testLengthBiasFail:                                                                                                                                                 

Annotations: There are 10 sets of 10 annotations of incremental size.

Segments: Segments are distributed uniformly over the complete workspace.

Workspace: Complete

Expected result: Longer annotations are called enriched when computing hypergeometric test
based on number of annotations that overlap.

Pass case: (default gat analysis)

Expected result: No annotations will be reported as enriched.

Chromosomal bias
================

Segments and annotations are not equally distributed amongst chromosomes, for example 
between autosomes and sex chromosomes. Ignoring the differences between chromosomes will
result in false depletion calls.

testChromosomalBiasFail:
* Annotations: There are 10 sets of 10 annotations of incremental size in 
  two workspaces. One workspace contains another set of 10 annotations. The workspaces are 
  then concatenated.
* Segments: Segments are distributed uniformly over the first workspace only.
* Workspace: Complete. 
* Expected result: The second set of annotations will be called depleted.

testChromosomalBiasPass:
* Workspace: Two workspaces
* Expected result: No annotations will be reported as depleted.

Composition bias
================

Segments and annotations might be associated with particular genomic features, 
such as G+C content. To correct for these, ischores can be defined.

Pass case:

The workspace is divided conceptually into 10 regions of similar size. 
The 10 regions are non-contiguous and are interleaved.

Annotations: There are 10 annotation sets of equal size that are distributed in the isochores 
such that there is always one annotation that is the prominent one (twice as large as all the others
combined).

Segments: Segments are distributed uniformly in workspace

Workspace: Complete workspace

Isochores: deactivated

Expected result: No enrichment is reported. Observed and sampled
segments will fall with equal proportion into the annotations.

Pass case:

Isochores: activated

Expected result: No enrichment is reported. Observed and sampled
segments will fall with equal proportion into the annotations.

Fail case:

Segments: Segments are distributed with different densities in
isochores.

Isochores: deactivated

Expected result: Some annotations will be reported enriched.

Pass case:

Isochores: activated

Expected result: No annotations will be reported enriched. The
different densities are accounted for by GAT.


NB:
* Simulations must be aware of edge effects. Intervals can't be placed
too regular.

* Simulations must test different densities. Especially I need to test
  if isochores and segments are of similar size or segments larger
  than isochores.

'''

import collections
import unittest
import gat
import numpy
import scipy.stats


def getRegularAnnotations(sizes=(100, 200, 400, 800,
                                 1600, 3200, 6400, 128000,
                                 256000, 512000),
                          nsegments = 20,
                          distance = 100):
    '''return annotations of size *size*. Each annotation has
    *nsegments* and segments are *distance* bases apart.

    Returns a dictionary of segment lists and start/end of 
    the annotations.
    '''

    annos = [[] for x in range(len(sizes))]

    start = distance
    for x in range(nsegments):
        for s, size in enumerate(sizes):
            annos[s].append((start, start + size))
            start += size + distance

    annotations = {}
    for x, size in enumerate(sizes):
        annotations["size%06i" % size] = gat.SegmentList(iter=annos[x],
                                                         normalize=True)
    return annotations, distance, start - distance


def getRegularSegments(workspace_size, size, density):
    '''return regular segments of size *size* and *density* in
    workspace.
    '''

    nsegments = int(workspace_size * density / size)
    rest = workspace_size * (1.0 - density)
    distance = rest // (nsegments + 2)
    assert distance > 0
    start = int(distance)
    s = []
    for x in range(nsegments):
        s.append((start, start + size))
        start += distance + size

    return gat.SegmentList(iter=s, normalize=True)


def getIncrementallySpacedSegments(workspace_size, size, density):
    '''return segments of size *size* and *density* in
    workspace. The inter-segment distance increases over the range.

    '''

    nsegments = int(workspace_size * density / size)
    rest = workspace_size * (1.0 - density)
    ngaps = nsegments + 2

    # sum of distances is total with increment.
    increment = rest // sum(range(ngaps))

    assert increment > 0

    distance = increment
    start = int(distance)
    s = []
    for x in range(nsegments):
        s.append((start, start + size))
        distance += increment
        start += distance + size

    return gat.SegmentList(iter=s, normalize=True)


def getRegularSegmentsInAnnotations(annotations, segment_size, segment_density):
    '''get a collection of segments overlapping with annotations.

    Segments are uniformly distributed over annotations. Longer annotations
    will have more segments. 

    Segments straddling annotations are truncated.
    '''

    merged = gat.SegmentList()
    for x, i in annotations.items():
        merged.extend(i)
    merged.normalize()

    workspace_size = merged.sum()
    nsegments = int(workspace_size * segment_density / segment_size)
    rest = workspace_size * (1.0 - segment_density)
    distance = rest // (nsegments + 2)

    # get regularly placed segments in a virtual workspace
    segments = getRegularSegments(
        workspace_size, segment_size, segment_density)

    # size of annotations
    n = 0
    m = 0
    s = []
    # negative for inter-segment segments
    overhang = distance
    annotation_start, annotation_end = merged[m]
    lannotation = annotation_end - annotation_start

    for start, end in segments:

        lsegment = end - start

        # place intergap segment, split if necessary
        while overhang > 0:

            if overhang < lannotation:
                annotation_start += overhang
                overhang = 0
                break

            overhang -= lannotation
            m += 1
            annotation_start, annotation_end = merged[m]
            lannotation = annotation_end - annotation_start

        lannotation = annotation_end - annotation_start

        # place segment - split if necessary
        while lsegment > 0:

            if lsegment < lannotation:
                s.append((annotation_start, annotation_start + lsegment))
                annotation_start += lsegment
                break

            s.append((annotation_start, annotation_end))

            lsegment -= lannotation
            m += 1
            annotation_start, annotation_end = merged[m]
            lannotation = annotation_end - annotation_start

        lannotation = annotation_end - annotation_start

    segments = gat.SegmentList(iter=s,
                               normalize=True)

    noverlap = segments.overlapWithSegments(merged)
    return segments


class GatTest(unittest.TestCase):
    '''test relative sizes of segments versus annotations.

    If the annotations cover the workspace completely, the enrichment in overlap
    in one annotation must be offset by the depletion in another annotation.

    This test iterates over a large range of possible values varying:

    1. the segment sizes
    2. the size of the annotations
    3. the arrangement of the annotations (single, versus several consecutive ones)
    '''

    counter = "NucleotideOverlap"

    nsamples_sampler = 3

    nsamples = 10

    variance = 0.2

    # P-Value of t-test.
    test_threshold = 0.05

    def setUp(self):
        self.outfile = open("validation_test.data", "w")
        self.bedfile = open("validation_test.bed", "w")

    def getSamples(self,
                   segments,
                   annotations,
                   workspace,
                   nsamples):

        sampler = gat.SamplerAnnotator(bucket_size=1, nbuckets=100000)

        if self.counter == "NucleotideOverlap":
            counter = gat.CounterNucleotideOverlap()
        elif self.counter == "SegmentOverlap":
            counter = gat.CounterSegmentOverlap()

        result = collections.defaultdict(list)

        for x in range(self.nsamples_sampler):
            r = gat.run(segments,
                        annotations,
                        workspace,
                        sampler,
                        counter,
                        num_samples=nsamples)

            for anno, v in r["default"].items():
                result[anno].append(v)

        return result

    def dumpSegments(self, segments, annotations, workspace, track):
        '''output segments, annotations and workspaces into a bed-formatted file.'''

        segments.save(self.bedfile, prefix="segments_", color="255,0,0")

        annotations.save(self.bedfile, prefix="annotations_", color="0,255,0")

        self.bedfile.write("track name=workspace_%s color=0,0,255\n" % track)
        for contig, v in workspace.items():
            for start, end in v:
                self.bedfile.write("%s\t%i\t%i\n" % (contig, start, end))

    def checkSignificant(self, result, fold_is_different):
        '''check if observed fold changes are significantly different/not different
        from 1.0.
        '''

        pvalues = []

        self.outfile.write("key\tmean\tpvalue\tmean_values\n")

        fail = []

        for key, rr in sorted(result.items()):

            folds = []
            for r in rr:
                pvalues.append(r.pvalue)
                folds.append(r.fold)

            mean_fold = numpy.mean(folds)

            # test if significantly different from 1
            tstat, pvalue = scipy.stats.ttest_1samp(folds, 1.0)

            self.outfile.write("%s\t%f\t%f\t%s\n" % (key, mean_fold, pvalue,
                                                     ",".join(map(str, folds))))

            # If p-value higher, we can not reject the
            # the null hypothesis that the mean of fold changes is
            # 1.0
            # Test not strictly applicable in both ways.
            if not fold_is_different:
                if pvalue < self.test_threshold:
                    fail.append("%s: expected fold change not to be different from 1.0. mean_fold=%f, pvalue=%f, folds=%s" % (
                        key, mean_fold, pvalue, folds))
            else:
                if pvalue >= self.test_threshold:
                    fail.append("%s: expected fold change to be different from 1.0. mean_fold=%f, pvalue=%f, folds=%s" % (
                        key, mean_fold, pvalue, folds))

        # check if pvalues are uniformly distributed between 0 and 1
        # a certain proportion of tests may pass
        print(pvalues)

        if fail:
            for f in fail:
                self.outfile.write("# failed: %s\n" % f)
            self.assertEqual(len(fail), 0, "%i tests failed" % len(fail))

    def printResults(self, result, section):

        self.outfile.write("section: %s\n" % section)
        self.outfile.write("%s\n" % "\t".join(gat.AnnotatorResult.headers))

        pvalues = []
        folds = []
        keys = sorted(result.keys())
        for key in keys:
            rr = result[key]
            f = []
            for r in rr:
                self.outfile.write("%s\n" % r)
                pvalues.append(r.pvalue)
                f.append(r.fold)

            folds.append(f)

    def addSingleIsochore(self, segments, annotations, workspace):
        '''add a single isochore to all segment lists.'''

        ss = gat.IntervalCollection("segment")
        ss.add("default", "chr1", segments)

        aa = gat.IntervalCollection("annotation")
        for key, x in annotations.items():
            aa.add(key, "chr1", x)

        ww = {"chr1": workspace}

        return ss, aa, ww

    def addIsochores(self, segments, annotations, workspace):
        '''add a single isochore to all segment lists.'''

        assert len(segments) == len(annotations) == len(workspace)

        ss = gat.IntervalCollection("segment")
        aa = gat.IntervalCollection("annotation")

        ww = {}

        chrom = 1
        for s, a, w in zip(segments, annotations, workspace):
            c = "chr%i" % chrom
            if s:
                ss.add("default", c, s)

            for key, x in a.items():
                if x:
                    aa.add(key, c, x)

            if w:
                ww[c] = w
            chrom += 1

        return ss, aa, ww

    def check(self,
              segments, annotations, workspace,
              label, fold_is_different):

        self.dumpSegments(segments, annotations, workspace, label)

        result = self.getSamples(
            segments, annotations, workspace, self.nsamples)

        self.printResults(result, label)
        self.checkSignificant(result, fold_is_different=False)


class ValidationTest(GatTest):

    segment_density = 0.1

    segment_size = 10

    def testDefault(self):

        annotations, start, end = getRegularAnnotations()

        workspacelist = gat.SegmentList(iter=[(start - 100, end + 100)],
                                        normalize=True)

        segmentlist = getRegularSegments(workspacelist.sum(),
                                         self.segment_size,
                                         self.segment_density)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testDefaultPass", fold_is_different=False)

    def testAscertainmentBiasFail1(self):

        segment_density = 0.1
        segment_size = 100

        annotations, start, end = getRegularAnnotations()

        # workspace twice as large as needed
        workspacelist = gat.SegmentList(iter=[(start - 100, 2 * end + 100)],
                                        normalize=True)

        segmentlist = getRegularSegments(workspacelist.sum(),
                                         self.segment_size,
                                         self.segment_density)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testAscertainmentBiasFail1",
                   fold_is_different=True)

    def testAscertainmentBiasFail2(self):

        annotations, start, end = getRegularAnnotations()

        workspacelist = gat.SegmentList(iter=[(start - 100, end + 100)],
                                        normalize=True)

        segmentlist = getRegularSegmentsInAnnotations(annotations,
                                                      self.segment_size,
                                                      self.segment_density)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testAscertainmentBiasFail2",
                   fold_is_different=True)

    def testAscertainmentBiasPass(self):

        annotations, start, end = getRegularAnnotations()

        workspacelist = gat.SegmentList(iter=[(start - 100, end + 100)],
                                        normalize=True)

        segmentlist = getRegularSegmentsInAnnotations(annotations,
                                                      self.segment_size,
                                                      self.segment_density)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testAscertainmentBiasPass",
                   fold_is_different=False)

    def testLengthBiasFail(self):

        annotations, start, end = getRegularAnnotations()

        workspacelist = gat.SegmentList(iter=[(start - 100, end + 100)],
                                        normalize=True)

        segmentlist = getRegularSegments(workspacelist.sum(),
                                         self.segment_size,
                                         self.segment_density)

    def testChromosomalBiasFail(self):

        annotations1, start, end = getRegularAnnotations()
        annotations2, start, end = getRegularAnnotations()

        # concatenate annotations
        annotations2.shift(workspacelist.sum())
        annotations1.extend(annotations2)

        workspacelist = gat.SegmentList(iter=[(start - 100, 2 * end + 100)],
                                        normalize=True)

        # segments only in first part of workspace
        segmentlist = getRegularSegments(end - start + 200,
                                         self.segment_size,
                                         self.segment_density)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testChromosomalBiasFail",
                   fold_is_different=True)

    def testChromosomalBiasPass(self):

        annotations1, start, end = getRegularAnnotations()
        annotations2, start, end = getRegularAnnotations()

        workspacelist = gat.SegmentList(iter=[(start - 100, end + 100)],
                                        normalize=True)

        # segments only in first part of workspace
        segmentlist = getRegularSegments(end - start + 200,
                                         self.segment_size,
                                         self.segment_density)

        ss, aa, ww = self.addIsochores((segmentlist, None),
                                       (annotations1, annotations2),
                                       (workspacelist, workspacelist))

        self.check(ss, aa, ww, "testChromosomalBiasPass",
                   fold_is_different=False)

    def testCompositionBiasFail(self):

        isochores = gat.SegmentList()
        for x in range(self.nisochores):

            start = x * self.isochore_size * self.nisochores

            i = start
            for x in range(self.nisochores):
                isochores.add(i, i + self.isochore_size)
                i += self.isochore_size

            # segments have increasing spacing within isochores
            s = getIncrementallySpacedSegments(self.nisochores * self.isochore_size,
                                               self.segment_size,
                                               self.segment_density)
            s.shift(start)
            segmentlist.extend(s)

            # annotations =


class SanityTest(GatTest):
    '''perform sanity checks.'''

    def testOverextendingIsochoresPass(self):
        '''
        annotations: covering full workspace
        isochores: covering full workspace
        workspace: covering full isochores 
        segments: uniformly distributed in isochores
        '''

        annotations = {}

        ws_size = 1000000
        nworkspaces = 20

        for x in range(nworkspaces):
            annotations["chr%i" % x] = gat.SegmentList(
                iter=[(x * ws_size, (x + 1) * ws_size)],
                normalize=True)

        workspacelist = gat.SegmentList(iter=[(0, nworkspaces * ws_size)],
                                        normalize=True)

        segmentlist = getRegularSegments(nworkspaces * ws_size, 1, 0.001)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testChromosomalBiasFail",
                   fold_is_different=False)

    def testOverextendingIsochoresFail(self):
        '''
        annotations: covering full workspace
        isochores: covering full workspace
        workspace: covering isochores partially
        segments: uniformly distributed in isochores
        '''

        annotations = {}

        ws_size = 1000000
        nworkspaces = 20

        for x in range(nworkspaces):
            annotations["chr%i" % x] = gat.SegmentList(
                iter=[(x * ws_size, (x + 1) * ws_size)],
                normalize=True)

        workspacelist = gat.SegmentList(iter=[(0, nworkspaces * ws_size)],
                                        normalize=True)

        segmentlist = getRegularSegments(nworkspaces * ws_size, 1, 0.001)

        ss, aa, ww = self.addSingleIsochore(
            segmentlist, annotations, workspacelist)

        self.check(ss, aa, ww, "testChromosomalBiasFail",
                   fold_is_different=False)

if __name__ == '__main__':
    unittest.main()
