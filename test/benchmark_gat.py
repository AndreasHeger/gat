'''test gat correctness
'''

import unittest
import random
import tempfile
import shutil
import os
import re
import gzip
import sys
import subprocess
import collections
import numpy
import math
import matplotlib.pyplot as plt

from GatSegmentList import SegmentList
import GatEngine
import gat

ANNOTATOR_DIR = "../annotator"
ANNOTATOR_CMD = '''java -Xmx8000M -cp %s/lib/commons-cli-1.0.jar:%s/lib/Annotator.jar app.Annotator -verbose 4''' % (
    ANNOTATOR_DIR, ANNOTATOR_DIR)


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also: 

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = numpy.r_[2 * x[0] - x[window_len:1:-1],
                 x, 2 * x[-1] - x[-1:-window_len:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='same')
    return y[window_len - 1:-window_len + 1]


def getDiff(a, b):
    if a == b == 0:
        return 0
    if a == 0:
        return 1
    return abs(float(a - b) / a)


def getPlotFilename(s, section=None):
    if section:
        filename = "%s_%s.png" % (re.sub("[ ()]", "", s), section)
    else:
        filename = "%s.png" % re.sub("[ ()]", "", s)
    filename = re.sub("__main__", "", filename)
    return os.path.join('benchmark', filename)


def getFilename(s):
    '''return filename for benchmark data.

    Benchmark data is stored in the _benchmark_ directory.
    '''
    filename = "%s" % re.sub("[ ()]", "", s)
    filename = re.sub("__main__", "", filename)
    return os.path.join('benchmark', filename)


class GatTest(unittest.TestCase):

    def shortDescription(self):
        return None


def writeAnnotatorSegments(outfile, segmentlist, section):

    if section == "workspace":
        prefix = "##Work"
    elif section == "segments":
        prefix = "##Seg"
    elif section == "annotations":
        prefix = "##Id\t0"

    contig = "chr1"
    os.write(outfile,
             "%s\t%s\t%s\n" % (
                 prefix,
                 contig,
                 "\t".join(["(%i,%i)" % (x, y) for x, y in segmentlist])))

    if section == "annotations":
        os.write(outfile, "##Ann\ttest\t0")

    os.close(outfile)


def writeAnnotatorSegments2(outfile, segmentlists, section):

    if section == "workspace":
        prefix = "##Work"
    elif section == "segments":
        prefix = "##Seg"
    elif section == "annotations":
        prefix = "##Id\t0"

    for contig, segmentlist in list(segmentlists.items()):
        os.write(outfile,
                 "%s\t%s\t%s\n" % (
                     prefix,
                     contig,
                     "\t".join(["(%i,%i)" % (x, y) for x, y in segmentlist])))

        if section == "annotations":
            os.write(outfile, "##Ann\ttest\t0")

    os.close(outfile)


def writeAnnotatorAnnotations(outfile, annotations):

    prefix = "##Id\t"
    nid = 0
    for track, segmentlists in annotations.items():
        nids = []
        for contig, segmentlist in list(segmentlists.items()):
            os.write(outfile,
                     "%s\t%i\t%s\t%s\n" % (
                         prefix,
                         nid,
                         contig,
                         "\t".join(["(%i,%i)" % x for x in segmentlist])))
            nids.append(nid)
            nid += 1

        os.write(outfile, "##Ann\t%s\t%s\n" %
                 (track, "\t".join(map(str, nids))))

    os.close(outfile)


def getAnnotatorSamples(segments,
                        annotations,
                        workspace,
                        sample_size):
    '''obtain results from annotator.'''
    fsegments, nsegments = tempfile.mkstemp()
    fworkspace, nworkspace = tempfile.mkstemp()
    fannotations, nannotations = tempfile.mkstemp()

    writeAnnotatorSegments2(fworkspace, workspace, "workspace")
    writeAnnotatorSegments2(fsegments, segments["default"], "segments")
    writeAnnotatorAnnotations(fannotations, annotations)

    statement = " ".join((ANNOTATOR_CMD,
                          "-seed %i" % random.getrandbits(16),
                          "-dumpsegments",
                          "-iterations %i" % sample_size,
                          # increment to avoid removing results with P=1.0
                          "-baseline 1.01",
                          "-annotation %s" % nannotations,
                          "-segments %s" % nsegments,
                          "-workspace %s" % nworkspace))

    process = subprocess.Popen(statement,
                               cwd=os.getcwd(),
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    # process.stdin.close()
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError("Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" %
                      (-process.returncode, stderr, statement))

    annotator_results = []

    intersection = SegmentList(clone=segments["default"]["chr1"])
    intersection.intersect(workspace["chr1"])
    segment_size = intersection.sum()

    # intersection = SegmentList( clone = annotations["default"]["chr1"] )
    # intersection.intersect( workspace["chr1"] )
    # annotation_size = intersection.sum()
    # print "segment_size", segment_size, annotation_size

    keep = False
    for line in stdout.split("\n"):

        # print line
        if line.startswith("#"):
            continue
        if line.startswith("Z"):
            keep = True
            continue
        if not keep:
            continue

        data = line.split("\t")
        if len(data) != 9:
            break
        z, fold, pvalue, observed, expected, lower95, upper95, stddev, annotation = data
        pvalue, observed, expected, lower95, upper95, stddev = \
            list(map(float, (pvalue, observed, expected, lower95, upper95, stddev)))

        try:
            fold = observed / expected
        except ZeroDivisionError:
            fold = 0.0

        r = GatEngine.DummyAnnotatorResult()
        r.track = "default"
        r.observed = observed * segment_size
        r.expected = expected * segment_size
        r.annotation = annotation[1:]
        r.fold = fold
        r.lower95 = lower95 * segment_size
        r.upper95 = upper95 * segment_size
        r.pvalue = pvalue
        r.qvalue = 1.0
        r.stddev = stddev * segment_size

        annotator_results.append(r)

    os.unlink(nsegments)
    os.unlink(nworkspace)
    os.unlink(nannotations)

    return annotator_results

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestSamplerLength(GatTest):
    '''sample from length distribution.'''
    ntests = 1000
    nsegments = 10000
    sampler = None

    def setUp(self):
        pass

    def testSamplingNormalDistribution(self):

        if not self.sampler:
            return

        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        self.segments = SegmentList(iter=[(x, x + numpy.random.randn() * 10.0 + 100.0)
                                          for x in range(0, 1000 * self.nsegments, 1000)],
                                    normalize=True)

        self.histogram = self.segments.getLengthDistribution(
            1, 1000 * self.nsegments)

        hs = self.sampler(self.histogram, 1)

        samples = [hs.sample() for x in range(0, 1 * self.nsegments)]

        delta = 1.0

        self.assertTrue(abs(numpy.mean(samples) - 100.0) < delta)
        self.assertTrue(abs(numpy.std(samples) - 10.0) < delta)

    def testSamplingSNPs(self):

        if not self.sampler:
            return

        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        self.segments = SegmentList(iter=[(x, x + 1)
                                          for x in range(0, 1000 * self.nsegments, 1000)],
                                    normalize=True)

        self.histogram = self.segments.getLengthDistribution(
            1, 1000 * self.nsegments)

        hs = self.sampler(self.histogram, 1)

        samples = [hs.sample() for x in range(0, 1 * self.nsegments)]

        self.assertAlmostEqual(numpy.mean(samples),
                               1.0,
                               places=0)

        self.assertAlmostEqual(numpy.std(samples),
                               0.0,
                               places=0)


class TestSamplerLengthFast(TestSamplerLength):

    sampler = GatEngine.HistogramSampler


class TestSamplerLengthSlow(TestSamplerLength):

    sampler = GatEngine.HistogramSamplerSlow

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestSegmentTrimming(GatTest):
    '''test for any biases in the trimming procedure.'''

    ntests = 10000

    def trim(self, segment_list, amount):
        '''trim amount from segment list.'''

        temp_sampler = GatEngine.SegmentListSampler(segment_list)

        # sample a position from current list of sampled segments
        start, end, overlap = temp_sampler.sample(1)

        segment_list.trim_ends(start,
                               amount,
                               numpy.random.randint(0, 2))

        return start

    def getSamples(self, segment_list, amount):

        samples = []
        starts = []
        for x in range(self.ntests):
            s = SegmentList(clone=segment_list)
            start = self.trim(s, amount)
            samples.append(s)
            starts.append(start)

        return samples, starts

    def getWorkspaceCounts(self,
                           segments,
                           amount):

        workspace = SegmentList(
            iter=[(segments.min(), segments.max()), ], normalize=True)

        nsegments = len(segments)
        workspace_size = workspace.sum()

        # segment_density = sample_length / workspace_size
        sample_length = segments.sum()

        expected = self.ntests * sample_length / workspace_size
        segment_density = sample_length / workspace_size

        samples, starts = self.getSamples(segments, amount)

        # print numpy.histogram( starts, bins = xrange( 0, workspace_size+1) )

        counts_within_workspace = getWorkspaceCounts(workspace,
                                                     segments,
                                                     samples,
                                                     filename=getPlotFilename(
                                                         str(self)),
                                                     expected=expected,
                                                     density=segment_density)

        d = abs(counts_within_workspace.mean() - expected) / float(expected)
        self.assertTrue(d < 0.1, "expected counts (%f) != sampled counts (%f" % (expected,
                                                                              counts_within_workspace.mean()))

        d = numpy.std(counts_within_workspace)

        return counts_within_workspace

    def testUniformDistance(self):
        '''test trimming with uniformly spaced segments
        '''

        segment_size = 10
        space = 10
        nsegments = 3
        amount = 4

        segments = SegmentList(iter=[(x, x + segment_size) for x in
                                     range(0,
                                            nsegments *
                                            (segment_size + space),
                                            segment_size + space)],
                               normalize=True)

        counts_within_workspace = self.getWorkspaceCounts(segments,
                                                          amount)

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestPositionSampling(GatTest):
    '''test segment position sampling (fixed length).'''

    ntests = 10000

    def getSamples(self, workspace, sampler, sample_length):

        samples = []
        for x in range(self.ntests):
            start, end, overlap = sampler.sample(sample_length)
            samples.append([(start, end)])

        return samples

    def getWorkspaceCounts(self,
                           workspace,
                           sampler,
                           sample_length):

        # +1 sample_length overhanging ends
        nsegments = len(workspace)
        workspace_size = workspace.sum() + (nsegments * (sample_length - 1))

        # segment_density = sample_length / workspace_size

        expected = self.ntests * sample_length / workspace_size
        segment_density = sample_length / workspace_size

        samples = self.getSamples(workspace, sampler, sample_length)

        counts_within_workspace = getWorkspaceCounts(workspace,
                                                     workspace,
                                                     samples,
                                                     filename=getPlotFilename(
                                                         str(self)),
                                                     expected=expected,
                                                     density=segment_density)

        d = abs(counts_within_workspace.mean() - expected) / float(expected)
        self.assertTrue(d < 0.1, "expected counts (%f) != sampled counts (%f" % (expected,
                                                                              counts_within_workspace.mean()))

        d = numpy.std(counts_within_workspace)

        return counts_within_workspace

    def testSingleWorkspace(self):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace_size = 10000
        sample_length = 4

        workspace = SegmentList(iter=[(0, workspace_size)],
                                normalize=True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(workspace,
                                                          sampler,
                                                          sample_length)

    def testTinyWorkspace(self):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace_size = 12
        sample_length = 4

        workspace = SegmentList(iter=[(0, workspace_size)],
                                normalize=True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(
            workspace, sampler, sample_length)

    def testSplitWorkspace(self):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace = SegmentList(iter=[(0, 1000), (9000, 10000)],
                                normalize=True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(
            workspace, sampler, 10)

    def testSplitWorkspace2(self):
        '''
        test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        # note if the the space is too small, there are problems
        # in counting the expected overlap.
        workspace = SegmentList(iter=[(0, 100), (200, 300)],
                                normalize=True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(workspace,
                                                          sampler,
                                                          2)

    def testMultipleWorkspaces(self):
        '''
        test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace = SegmentList(iter=((x, x + 100) for x in range(0, 10000, 1000)),
                                normalize = True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(
            workspace, sampler, 10)

    def testSNPPositionSampling(self):
        '''
        test if we sample the exactly right amount of nucleoutides.'''

        workspace = SegmentList(iter=((x, x + 100) for x in range(0, 10000, 1000)),
                                normalize = True)

        sampler = GatEngine.SegmentListSampler(workspace)

        counts_within_workspace = self.getWorkspaceCounts(
            workspace, sampler, 1)

##############################################################################
##############################################################################
##############################################################################
##############################################################################


def getWorkspaceCounts(workspace,
                       segments,
                       samples,
                       filename,
                       expected=None,
                       density=0):
    '''compute sample counts within workspace.

    returns array of size workspace with counts.

    Also outputs plots.
    '''

    l = workspace.max()
    counts = numpy.zeros(l, numpy.int)
    ntests = len(samples)

    segment_sizes = []
    nsegments = []
    starts = numpy.zeros(l + 1, numpy.int)
    ends = numpy.zeros(l + 1, numpy.int)

    for sample_id, s in enumerate(samples):

        if len(s) == 0:
            continue

        for start, end in s:
            start = max(0, start)
            end = min(end, l)
            counts[start:end] += 1
            starts[start] += 1
            ends[end] += 1

        ss = [x[1] - x[0] for x in s]

        segment_sizes.append((numpy.std(ss),
                              min(ss),
                              max(ss),
                              numpy.mean(ss),
                              numpy.median(ss)))
        nsegments.append(len(ss))

    counts_within_workspace = []
    for start, end in workspace:
        counts_within_workspace.extend(counts[start:end])

    if l > 10:
        dx = 10
    else:
        dx = 1

    counts_within_workspace = numpy.array(
        counts_within_workspace, dtype=numpy.int)
    newy = smooth(counts_within_workspace, window_len=dx)

    plt.figure(figsize=(14, 8), dpi=80)
    plt.subplots_adjust(right=0.7)
    # plt.axes( [0.1,0.1,0.51,0.5] )

    plt.subplot("511")
    plt.plot(range(len(counts_within_workspace)),
             counts_within_workspace,
             '.',
             label="coverage")

    plt.plot(range(len(counts_within_workspace)), newy, '-',
             label="smooth (%i)" % dx)

    plt.title("%s : density = %6.4f" % (filename, density))
    # : nworkspaces=%i, sample_length=%i" % ( filename, len(workspace), len(samples) ) )
    plt.xlabel("position")
    plt.ylabel("counts")
    if expected:
        d = expected * 0.1
        plt.plot(range(len(counts_within_workspace)), [expected] * len(counts_within_workspace),
                 '-', label="expected")
        plt.plot(range(len(counts_within_workspace)), [expected - d] * len(counts_within_workspace),
                 '--')
        plt.plot(range(len(counts_within_workspace)), [expected + d] * len(counts_within_workspace),
                 '--')
    plt.legend(loc=(1.03, 0.2))

    plt.subplot("512")
    segment_sizes.sort()
    segment_sizes = list(zip(*segment_sizes))
    plt.plot(segment_sizes[0], label="stddev")
    plt.plot(segment_sizes[1], label="min")
    plt.plot(segment_sizes[2], label="max")
    plt.plot(segment_sizes[3], label="mean")
    plt.plot(segment_sizes[4], label="median")
    plt.legend(loc=(1.03, 0))
    plt.xlabel("sample")
    plt.ylabel("segment size")

    plt.subplot("513")
    plt.plot(starts, label="starts")
    plt.plot(ends, label="ends")
    plt.legend(loc=(1.03, 0.2))
    plt.xlabel("position")
    plt.ylabel("counts")

    plt.subplot("514")
    nsegments.sort()
    plt.plot(nsegments, label="nsegments")
    plt.legend(loc=(1.03, 0.2))
    plt.xlabel("sample")
    plt.ylabel("counts")
    plt.axhline(len(segments), color='r', ls='--')

    # ouptput test configuration
    plt.subplot("515")
    mi, ma = workspace.min(), workspace.max()

    extent = float(ma - mi)
    for start, end in segments:
        plt.axhspan(0.28, 0.3, start / extent, end / extent)

    for start, end in workspace:
        plt.axhspan(0.4, 0.42, start / extent, end / extent, color="r")

    plt.ylim(0.1, 0.5)
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.savefig(filename)

    return counts_within_workspace

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestSegmentSamplingSamplerAnnotator(GatTest):

    ntests = 1000

    # is as expected. Set to False for samplers
    # which do not return the exact number of nucleotides
    check_nucleotides = True

    # check average coverage, set to false for samples
    check_average_coverage = True
    # check for uniform covage
    check_uniform_coverage = True

    # check segment size distribution
    # turned off as segment length will v ary
    # in testFullWorkspace
    check_length = False

    def setUp(self):
        self.sampler = GatEngine.SamplerAnnotator()

    def getSample(self, segments, workspace):
        '''return a single sample'''
        return self.sampler.sample(segments, workspace)

    def getSamples(self, segments, workspace):
        '''return a list of samples.'''
        samples = []

        for x in range(self.ntests):
            sample = self.getSample(segments, workspace)
            samples.append(sample)

        return samples

    def checkLength(self, samples, segments, workspace):
        '''check length distribution of segments.'''

        # check segment lengths
        l = [x[1] - x[0] for x in segments]
        values_input = min(l), max(l), numpy.mean(l), numpy.std(l)

        for i, sample in enumerate(samples):
            l = [x[1] - x[0] for x in sample]
            values_sample = min(l), max(l), numpy.mean(l), numpy.std(l)

            for val, inp, samp in zip(("min", "max", "mean", "std"),
                                      values_input,
                                      values_sample):
                d = abs(inp - samp) / float(inp)
                self.assertTrue(d < 0.1,
                             "segment length distribution in sample %i: expected %s (%f) != observed %s (%f)" %
                             (i, val, inp, val, samp))

    def checkSanity(self, samples, segments, workspace):
        '''sanity check on segments.'''

        for i, sample in enumerate(samples):
            s = SegmentList(clone=sample)
            s.filter(workspace)
            self.assertEqual(len(s), len(
                sample), "sample %i contains segments outside of workspace %s: %s" % (i, workspace, sample))

    def checkSample(self, samples, segments, workspace):
        '''check if sample corresponds to expectation.

        segment_length: total length of segments (nsegments * segment_length)
        '''

        filename = getPlotFilename(str(self))

        self.assertEqual(self.ntests, len(samples))

        # compute expected coverage
        # check if coverage is uniform over workspace
        # expected coverage
        # number_of_samples * segment_length / workspace_length
        # modifier:
        # can be simplified
        nsegments = float(len(segments))
        tmp = SegmentList(clone=workspace)
        tmp.intersect(segments)
        # number of nucleotides in segments overlapping the workspace
        segment_overlap = tmp.sum()
        # density of segment nucleotides in workspace
        segment_density = float(segment_overlap) / workspace.sum()

        segment_length = segment_overlap / nsegments
        expected = len(samples) * segment_overlap / float(workspace.sum())
        #   float(workspace.sum()) / (workspace.sum() + segments.sum()  * len(workspace) )

        # compute actual coverage counts and plot
        counts_within_workspace = getWorkspaceCounts(workspace,
                                                     segments,
                                                     samples,
                                                     filename,
                                                     expected=expected,
                                                     density=segment_density)

        # for x in samples: print str(x)

        # check if each sample has the correct number of nucleotides
        sums = [x.sum() for x in samples]

        self.checkSanity(samples, segments, workspace)

        if self.check_length:
            self.checkLength(samples, segments, workspace)

        if self.check_nucleotides:
            for x, s in enumerate(samples):
                tmp = SegmentList(clone=workspace, normalize=True)
                s.normalize()
                tmp.intersect(s)
                ovl = tmp.sum()

                self.assertEqual(ovl, segment_overlap,
                                 "incorrect number of nucleotides in sample %i, got %i, expected %i\nsample=%s" %
                                 (x, ovl, segment_overlap, samples[x]))

        if self.check_average_coverage:
            # check average coverage
            d = abs(counts_within_workspace.mean() -
                    expected) / float(expected)

            self.assertTrue(d < 0.1, "expected counts (%f) != sampled counts (%f)" % (expected,
                                                                                   counts_within_workspace.mean()))

        if self.check_uniform_coverage:
            # check for uniform coverage
            stddev = numpy.std(counts_within_workspace)
            d = stddev / float(expected)

            self.assertTrue(d < 0.1, "coverage variation too large : stddev (%f) / %f = %f > 0.01" %
                         (stddev,
                          expected,
                          d))

    def testSegmentedWorkspaceSmallGap(self):
        '''
        test sampling within a segmented workspace.

        The gap between adjacent workspaces is 10 nucleotides
        to avoid the merging of intervals. The gap is thus
        much smaller than the segment size.
        '''
        nsegments = 10
        segment_size = 100

        workspace = SegmentList(iter=((x, x + 990)
                                      for x in range(0, 1000 * nsegments, 1000)),
                                normalize = True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, 1000 * nsegments, 1000)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspacePartiallyOverlappingSegments(self):
        '''
        test sampling within a segmented workspace where
        segments overlap only partially.

        w:   ----    ----    ----
        s: ----    ----    ----

        '''
        nsegments = 10
        segment_size = 100

        half = segment_size // 2
        workspace = SegmentList(iter=((x, x + segment_size)
                                      for x in range(half,
                                                     nsegments * segment_size,
                                                     2 * segment_size)),
                                normalize = True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0,
                                                    nsegments * 2 * segment_size,
                                                    2 * segment_size)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspaceSmallGapUnequalSides(self):
        '''
        test sampling within a segmented workspace.

        The workspace has two segments. The gap between
        segments is half a segment size.

        The gap between adjacent workspaces is 5 nucleotides.
        '''
        nsegments = 1
        segment_size = 50

        workspace = SegmentList(iter=((0, 1 * segment_size),
                                      (int(1.5 * segment_size), 2 * segment_size)),
                                normalize = True)
        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, 10 * nsegments, 20)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspaceSmallGapEqualSides(self):
        '''
        test sampling within a segmented workspace.

        The workspace has two segments. The gap between
        segments is half a segment size.

        The gap between adjacent workspaces is 5 nucleotides.
        '''
        nsegments = 1
        segment_size = 50
        gap = 5
        workspace = SegmentList(iter=((0, segment_size),
                                      (segment_size + gap, 2 * segment_size + gap)),
                                normalize = True)
        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, 10 * nsegments, 20)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspaceSmallGapEqualSidesManySegments(self):
        '''
        test sampling within a segmented workspace.

        The workspace has two segments. The gap between
        segments is half a segment size.

        The gap between adjacent workspaces is 5 nucleotides.
        '''
        nsegments = 10
        segment_size = 5
        gap = 5
        half = nsegments * segment_size
        workspace = SegmentList(iter=((0, half),
                                      (half + gap, 2 * half + gap)),
                                normalize = True)
        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, half, 2 * segment_size)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspaceLargeGap(self):
        '''
        test sampling within a segmented workspace.

        The gap is the size of a segment.
        '''
        nsegments = 10
        segment_size = 100
        workspace = SegmentList(iter=((x, x + 900)
                                      for x in range(0, 1000 * nsegments, 1000)),
                                normalize = True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, 1000 * nsegments, 1000)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSingleWorkspace(self):
        '''
        test sampling within a single continuous workspace.'''
        nsegments = 10
        segment_size = 100
        workspace_size = 10000

        workspace = SegmentList(iter=[(0, workspace_size)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, 1000 * nsegments, 1000)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSingleWorkspaceWithOffset(self):
        '''
        test sampling within a single continuous workspace.'''
        nsegments = 10
        segment_size = 100
        workspace_size = 10000
        offset = 10000

        workspace = SegmentList(iter=[(offset, offset + workspace_size)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(offset, offset + 1000 * nsegments, 1000)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSingleWorkspaceSingleSegment(self):
        '''
        test sampling within a single continuous workspace.'''
        nsegments = 1
        segment_size = 1000
        workspace_size = 10000

        workspace = SegmentList(iter=[(0, workspace_size)],
                                normalize=True)

        segments = SegmentList(iter=[(workspace_size // 2 - segment_size // 2,
                                      workspace_size // 2 + segment_size // 2)],
                               normalize=True)

        samples = self.getSamples(segments, workspace)
        self.checkSample(samples, segments, workspace)

    def testFullWorkspace(self):
        '''
        test sampling within a single continuous workspace
        and a single segment.

        The segment is larger than the workspace
        '''
        nsegments = 1
        segment_size = 200
        offset = 1

        workspace = SegmentList(iter=[(0, 100)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, offset * nsegments, offset)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)
        self.checkSample(samples, segments, workspace)

    def testSmallWorkspace(self):
        '''
        test sampling within a single continuous workspace.
        Segments are half the size of the workspace
        '''
        nsegments = 1
        segment_size = 50
        offset = 1

        workspace = SegmentList(iter=[(0, 100)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, offset * nsegments, offset)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testTinyWorkspace(self):
        '''
        test sampling within a single small continuous workspace.

        '''
        nsegments = 1
        segment_size = 4
        offset = 1

        workspace = SegmentList(iter=[(0, 12)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, offset * nsegments, offset)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)
        self.checkSample(samples, segments, workspace)

    def testSmallWorkspaceManySegments(self):
        '''
        test sampling within a single continuous workspace.
        Segments half the size of the workspace
        '''
        nsegments = 10
        segment_size = 5

        workspace = SegmentList(iter=[(0, 100)],
                                normalize=True)

        segments = SegmentList(iter=((x,  x + segment_size)
                                     for x in range(0, segment_size * nsegments * 2, 2 * segment_size)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    def testSegmentedWorkspace2x(self):
        '''
        test if workspace segments are only twice the size of segmments.
        '''

        nsegments = 10
        segment_size = 100

        workspace = SegmentList(iter=((x, x + 2 * segment_size)
                                      for x in range(0, 1000 * nsegments, 1000)),
                                normalize = True)
        segments = SegmentList(iter=((x, x + segment_size)
                                     for x in range(0, 1000 * nsegments, 1000)),
                               normalize = True)

        samples = self.getSamples(segments, workspace)

        self.checkSample(samples, segments, workspace)

    # def testVariableLengthSampling( self ):

    #     nsegments = 10000
    #     # create normaly distributed lengths of mean 100.0 and sigma = 10.0
    #     segments = SegmentList( iter = [ (x, x + numpy.random.randn() * 10.0 + 100.0  ) \
    #                                              for x in range(0, 1000 * nsegments, 1000) ],
    #                                 normalize = True )

    #     histogram = segments.getLengthDistribution( 1, 1000 * nsegments )
    #     hs = GatEngine.HistogramSampler( histogram, 1 )

    #     samples = [hs.sample() for x in range(0, 1 * nsegments )]

    #     self.assertAlmostEqual( numpy.mean(samples),
    #                             100.0,
    #                             places = 0 )

    #     self.assertAlmostEqual( numpy.std(samples),
    #                             10.0,
    #                             places = 0 )


class TestSegmentSamplingSamplerSegments(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = False
    check_average_coverage = False
    check_uniform_coverage = False

    def setUp(self):
        self.sampler = GatEngine.SamplerSegments()


class TestSegmentSamplingSamplerBruteForce(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = True
    check_average_coverage = True
    check_uniform_coverage = True

    def setUp(self):
        self.sampler = GatEngine.SamplerBruteForce()


class TestSegmentSamplingSamplerLocalPermutation(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = True
    check_average_coverage = True
    check_uniform_coverage = True

    def setUp(self):
        self.sampler = GatEngine.SamplerLocalPermutation()


class TestSegmentSamplingSamplerGlobalPermutation(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = True
    check_average_coverage = True
    check_uniform_coverage = True

    def setUp(self):
        self.sampler = GatEngine.SamplerGlobalPermutation()


class TestSegmentSamplingSamplerUniform(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = False
    check_average_coverage = False
    check_uniform_coverage = False

    def setUp(self):
        self.sampler = GatEngine.SamplerUniform(increment=1)


class TestSegmentSamplingSamplerShift(TestSegmentSamplingSamplerAnnotator):

    check_nucleotides = False
    check_average_coverage = False
    check_uniform_coverage = False

    def setUp(self):
        self.sampler = GatEngine.SamplerShift(radius=5)


class TestSegmentSamplingSamplerTheAnnotator(TestSegmentSamplingSamplerAnnotator):
    '''use annotator to sample segments.'''

    def getSamples(self,
                   segments,
                   workspace):

        fsegments, nsegments = tempfile.mkstemp()
        fworkspace, nworkspace = tempfile.mkstemp()
        fannotations, nannotations = tempfile.mkstemp()

        writeAnnotatorSegments(fsegments, segments, "segments")
        writeAnnotatorSegments(fworkspace, workspace, "workspace")
        writeAnnotatorSegments(fannotations, workspace, "annotations")

        statement = " ".join((ANNOTATOR_CMD,
                              "-seed %i" % random.getrandbits(16),
                              "-dumpsegments",
                              "-iterations %i" % self.ntests,
                              # increment to avoid removing results with P=1.0
                              "-baseline 1.01",
                              "-annotation %s" % nannotations,
                              "-segments %s" % nsegments,
                              "-workspace %s" % nworkspace))

        process = subprocess.Popen(statement,
                                   cwd=os.getcwd(),
                                   shell=True,
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # process.stdin.close()
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise OSError("Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" %
                          (-process.returncode, stderr, statement))

        samples = []
        for line in stdout.split("\n"):
            if not line.startswith("##Segments"):
                continue
            data = line.split("\t")[2:]
            coords = [list(map(int, x[1:-1].split(","))) for x in data]
            # print "annotator", coords
            samples.append(SegmentList(iter=coords))

        os.unlink(nsegments)
        os.unlink(nworkspace)
        os.unlink(nannotations)

        return samples

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestStatsSNPSampling(GatTest):
    '''test Stats by running a SNP (1-size interval) analysis.

    For SNPs, the hypergeometric distribution applies.
    '''

    sample_size = 10

    def check(self, workspace, annotations, segments):

        workspace_size = workspace["chr1"].sum()

        sampler = GatEngine.SamplerAnnotator(
            bucket_size=1, nbuckets=workspace_size)

        counters = [GatEngine.CounterNucleotideOverlap()]
        workspace_generator = GatEngine.UnconditionalWorkspace()

        # print segments["default"]["chr1"]
        # print workspace["chr1"]

        annotator_results = gat.run(segments,
                                    annotations,
                                    workspace,
                                    sampler,
                                    counters,
                                    workspace_generator=workspace_generator,
                                    num_samples=self.sample_size)

        outfile = sys.stdout

        self.assertEqual(workspace_size, workspace["chr1"].sum())

        tmp = SegmentList(clone=segments["default"]["chr1"])
        tmp.intersect(workspace["chr1"])
        segment_size = tmp.sum()

        anno_mean, anno_pvalue, anno_std = [], [], []
        dist_mean, dist_pvalue, dist_std = [], [], []

        annos = sorted(annotations.keys())

        for annotation in annos:

            r = [x for x in annotator_results if x.annotation == annotation]
            assert len(r) == 1
            result = r[0]

            # print workspace["chr1"]
            # print annotations[annotation]["chr1"]

            intersection = SegmentList(clone=annotations[annotation]["chr1"])
            intersection.intersect(workspace["chr1"])

            annotation_size = intersection.sum()

            # expected overlap for sampling with replacement
            expected_with = annotation_size / float(workspace_size)

            if annotation_size < workspace_size:

                # test sampling (# of expected)
                # sampling without replacement follows hypergeometric distribution
                # good = annotation_size
                # bad = workspace - annotation_size
                hyper = numpy.random.hypergeometric(annotation_size,
                                                    workspace_size - annotation_size,
                                                    segment_size,
                                                    self.sample_size)
                hyper.sort()

                # m = annotation_size # (white balls)
                # N = workspace_size # (total balls)
                # n = segment_size # (balls taken)
                # variance = float(n * m * ( N - n ) * ( N -m )) / (N * N * (N - 1 )  )

                # expected overlap for sampling without replacement
                expected_without = hyper.mean()
                error = hyper.std() * 4  # / segment_size

                expected_std = hyper.std()

                expected_pvalue = GatEngine.getNPTwoSidedPValue(
                    hyper, result.observed)

                # print "\t".join( map(str, (result,
                #                            expected_without,
                #                            expected_std,
                #                            expected_pvalue,
                #                            workspace_size,
                #                            annotation_size) ) )

                # for small sample size there might be no positive samples
                if error == 0 and segment_size < 3:
                    continue
            else:
                # annotations = workspace
                expected_without = segment_size
                expected_std = 0
                expected_error = 0
                expected_pvalue = 1.0
                error = 0.1

            self.assertTrue(abs(result.expected - expected_without) < error,
                         "simulated results deviate from hypergeometric expectation for annotation `%s`: sizes(seg=%i/anno=%i/work=%i) observed=%f, expected=%f (%f, margin=%f)" %
                         (annotation,
                          segment_size,
                          annotation_size,
                          workspace_size,
                          result.expected, expected_without,
                          result.expected - expected_without,
                          error))

            anno_mean.append(result.expected)
            anno_std.append(result.stddev)
            anno_pvalue.append(result.pvalue)

            dist_mean.append(expected_without)
            dist_std.append(expected_std)
            dist_pvalue.append(expected_pvalue)

            # plt.figure()
            # hist1, bins1 = numpy.histogram( hyper, new = True, bins=xrange(0, segment_size+10 ))
            # hist2, bins2 = numpy.histogram( result.samples, new = True, bins=bins1 )
            # plt.title( "%i - %i - %i" % (segment_size,
            #                              annotation_size,
            #                              workspace_size) )
            # plt.plot( bins1[:-1], hist1, label = "dist" )
            # plt.plot( bins2[:-1], hist2, label = "simulated" )
            # plt.legend()
            # plt.show()

        plt.figure()
        plt.subplot(221)
        plt.scatter(anno_mean, dist_mean)
        plt.xlabel("simulated - mean")
        plt.ylabel("distribution - mean")
        plt.plot(anno_mean, anno_mean, "b")
        plt.subplot(222)
        plt.scatter(anno_std, dist_std)
        plt.xlabel("simulated - std")
        plt.ylabel("distribution - std")
        plt.plot(anno_std, anno_std, "b")
        plt.subplot(223)
        plt.scatter(anno_pvalue, dist_pvalue)
        plt.xlabel("simulated - pvalue")
        plt.ylabel("distribution - pvalue")
        plt.plot(anno_pvalue, anno_pvalue, "b")
        plt.savefig(os.path.join('benchmark', "test_%s.png" %
                                 re.sub("[ ()]", "", str(self))))

    def testSingleSNP(self):

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspace_size = 1000
        nsnps = 1

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        segments.add("default", "chr1", SegmentList(iter=[(x, x + 1) for x in range(0, nsnps)],
                                                    normalize=True))

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(1, 100, 2):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(0, y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testMultipleSNPsFullOverlap(self):

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspace_size = 1000
        nsnps = 10

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # 10 snps
        segments.add("default", "chr1", SegmentList(iter=[(x, x + 1) for x in range(0, nsnps)],
                                                    normalize=True))

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(10, 110, 5):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(0, y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testMultipleSNPsPartialOverlap(self):
        '''test with multiple snps and decreasing
        amount of overlap with annotations.

        Tests if p-values are computed correctly.
        '''
        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspace_size = 1000
        nsnps = 100

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # 10 snps
        segments.add("default", "chr1", SegmentList(iter=[(x, x + 1) for x in range(0, nsnps)],
                                                    normalize=True))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, nsnps, 1):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(y, nsnps + y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testIntervalsPartialOverlap(self):
        '''test with intervals with 
        increasing amount of overlap.

        '''
        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspace_size = 1000

        size = 100

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # segment of size 10
        segments.add("default", "chr1", SegmentList(iter=[(0, size), ],
                                                    normalize=True))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, size):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(y, size + y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testWorkspaces(self):
        '''
        input:
            workspace = 500 segments of size 1000, separated by a gap of 1000
            annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
            segments = a SNP every 100 bp

        output:
        '''

        workspace_size = 1000000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1",
                       SegmentList(iter=[(x, x + 1000) for x in range(0, workspace_size, 2000)],
                                   normalize=True))

        # SNPs every 100bp
        segments.add("default", "chr1",
                     SegmentList(iter=[(x, x + 1) for x in range(0, workspace_size, 100)],
                                 normalize=True))

        start = 0
        annotations.add("%03i" % start, "chr1",
                        SegmentList(iter=[(0, workspace_size)],
                                    normalize=True))

        self.check(workspaces["default"], annotations, segments)

    def testFullAnnotation(self):

        workspace_size = 1000000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1",
                       SegmentList(iter=[(x, x + 1000) for x in range(0, workspace_size, 2000)],
                                   normalize=True))

        # SNPs every 100bp
        segments.add("default", "chr1",
                     SegmentList(iter=[(x, x + 1) for x in range(0, workspace_size, 100)],
                                 normalize=True))

        size = 1000
        for start in range(0, size, 100):
            annotations.add("%03i" % start, "chr1",
                            SegmentList(iter=[(start + x, start + x + size)
                                              for x in range(0, workspace_size, 2000)],
                                        normalize=True))

        self.check(workspaces["default"], annotations, segments)

##############################################################################
##############################################################################
##############################################################################
##############################################################################


class TestStatsGat(GatTest):
    '''
    '''

    sample_size = 10

    def getSamples(self,
                   segments,
                   annotations,
                   workspace):

        workspace_size = workspace["chr1"].sum()

        sampler = GatEngine.SamplerAnnotator(
            bucket_size=1, nbuckets=workspace_size)

        counters = [GatEngine.CounterNucleotideOverlap()]

        workspace_generator = GatEngine.UnconditionalWorkspace()

        return gat.run(segments,
                       annotations,
                       workspace,
                       sampler,
                       counters,
                       workspace_generator=workspace_generator,
                       num_samples=self.sample_size)

    def check(self, workspace, annotations, segments):

        workspace_size = workspace["chr1"].sum()

        # print segments["default"]["chr1"]
        # print workspace["chr1"]

        annotator_results = self.getSamples(segments,
                                            annotations,
                                            workspace)

        annotator_results = GatEngine.IOTools.flatten(annotator_results)

        outfile = sys.stdout

        self.assertEqual(workspace_size, workspace["chr1"].sum())

        tmp = SegmentList(clone=segments["default"]["chr1"])
        tmp.intersect(workspace["chr1"])
        segment_size = tmp.sum()

        anno_mean, anno_pvalue, anno_std = [], [], []
        dist_mean, dist_pvalue, dist_std = [], [], []

        annos = sorted(annotations.keys())

        for annotation in annos:
            r = [x for x in annotator_results if x.annotation == annotation]
            assert len(r) == 1
            result = r[0]

            intersection = SegmentList(clone=annotations[annotation]["chr1"])
            intersection.intersect(workspace["chr1"])

            annotation_size = intersection.sum()

            # expected overlap for sampling with replacement
            expected_with = annotation_size / float(workspace_size)

            if annotation_size < workspace_size:

                # test sampling (# of expected)
                # sampling without replacement follows hypergeometric distribution
                # good = annotation_size
                # bad = workspace - annotation_size
                hyper = numpy.random.hypergeometric(annotation_size,
                                                    workspace_size - annotation_size,
                                                    segment_size,
                                                    self.sample_size)
                hyper.sort()

                # m = annotation_size # (white balls)
                # N = workspace_size # (total balls)
                # n = segment_size # (balls taken)
                # variance = float(n * m * ( N - n ) * ( N -m )) / (N * N * (N - 1 )  )

                # expected overlap for sampling without replacement
                expected_without = hyper.mean()
                error = hyper.std() * 4  # / segment_size

                expected_std = hyper.std()

                expected_pvalue = GatEngine.getNPTwoSidedPValue(
                    hyper, result.observed)

                # print "\t".join( map(str, (result,
                #                            expected_without,
                #                            expected_std,
                #                            expected_pvalue,
                #                            workspace_size,
                #                            annotation_size) ) )

                # for small sample size there might be no positive samples
                if error == 0 and segment_size < 3:
                    continue
            else:
                # annotations = workspace
                expected_without = segment_size
                expected_std = 0
                expected_error = 0
                expected_pvalue = 1.0
                error = 0.1

            # print annotation_size, workspace_size - annotation_size,
            # segment_size, self.sample_size, result.observed, result.expected,
            # result.pvalue, expected_pvalue, result.stddev, expected_std

            self.assertTrue(abs(result.expected - expected_without) < error,
                         "simulated results deviates from hypergeometric expectation for annotation `%s`: sizes(seg=%i/anno=%i/work=%i) observed=%f, expected=%f (%f, margin=%f)" %
                         (annotation,
                          segment_size,
                          annotation_size,
                          workspace_size,
                          result.expected, expected_without,
                          result.expected - expected_without,
                          error))

            anno_mean.append(result.expected)
            anno_std.append(result.stddev)
            anno_pvalue.append(result.pvalue)

            dist_mean.append(expected_without)
            dist_std.append(expected_std)
            dist_pvalue.append(expected_pvalue)

            # plt.figure()
            # hist1, bins1 = numpy.histogram( hyper, new = True, bins=xrange(0, segment_size+10 ))
            # hist2, bins2 = numpy.histogram( result.samples, new = True, bins=bins1 )
            # plt.title( "%i - %i - %i" % (segment_size,
            #                              annotation_size,
            #                              workspace_size) )
            # plt.plot( bins1[:-1], hist1, label = "dist" )
            # plt.plot( bins2[:-1], hist2, label = "simulated" )
            # plt.legend()
            # plt.show()

        plt.figure()
        plt.subplot(221)
        plt.scatter(dist_mean, anno_mean)
        plt.ylabel("simulated - mean")
        plt.xlabel("distribution - mean")
        plt.plot(anno_mean, anno_mean, "b")
        plt.subplot(222)
        plt.scatter(dist_std, anno_std)
        plt.ylabel("simulated - std")
        plt.xlabel("distribution - std")
        plt.plot(anno_std, anno_std, "b")
        plt.subplot(223)
        plt.scatter(dist_pvalue, anno_pvalue)
        plt.ylabel("simulated - pvalue")
        plt.xlabel("distribution - pvalue")
        plt.plot(anno_pvalue, anno_pvalue, "b")
        plt.savefig(getPlotFilename(str(self)))

    def testSingleSNP(self):
        '''test stats using a single SNP

        The counts should follow a hypergeometric distribution.
        '''

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspace_size = 1000

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
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

        self.check(workspace, annotations, segments)

    def testMultipleSNPsFullOverlap(self):

        nsnps = 10
        workspace_size = 1000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # 10 snps
        segments.add("default", "chr1", SegmentList(iter=[(x, x + 1) for x in range(0, nsnps * 2, 2)],
                                                    normalize=True))

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(10, 110, 5):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(0, y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testMultipleSNPsPartialOverlap(self):
        '''test with multiple snps and decreasing
        amount of overlap with annotations.

        Tests if p-values are computed correctly.
        '''
        workspace_size = 1000
        nsnps = 100

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # 10 snps
        segments.add("default", "chr1", SegmentList(iter=[(x, x + 1) for x in range(0, nsnps * 2, nsnps)],
                                                    normalize=True))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, nsnps * 2, 2):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(y, nsnps + y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testIntervalsPartialOverlap(self):
        '''
        test with intervals with increasing amount of overlap.
        '''

        workspace_size = 1000
        segment_size = 100

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1", SegmentList(iter=[(0, workspace_size), ],
                                                      normalize=True))
        workspace = workspaces["default"]

        # segment of size segment_size
        segments.add("default", "chr1", SegmentList(iter=[(0, segment_size), ],
                                                    normalize=True))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, segment_size):
            annotations.add("%03i" % y, "chr1",
                            SegmentList(iter=[(y, segment_size + y), ],
                                        normalize=True))

        self.check(workspace, annotations, segments)

    def testWorkspaces(self):
        '''
        input:
            workspace = 500 segments of size 1000, separated by a gap of 1000
            annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
            segments = a SNP every 100 bp

        output:
        '''

        workspace_size = 1000000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1",
                       SegmentList(iter=[(x, x + 1000) for x in range(0, workspace_size, 2000)],
                                   normalize=True))

        # SNPs every 100bp
        segments.add("default", "chr1",
                     SegmentList(iter=[(x, x + 1) for x in range(0, workspace_size, 100)],
                                 normalize=True))

        size = 1000
        for start in range(0, size, 100):
            annotations.add("%03i" % start, "chr1",
                            SegmentList(iter=[(start + x, start + x + size)
                                              for x in range(0, workspace_size, 2000)],
                                        normalize=True))

        self.check(workspaces["default"], annotations, segments)

    def testCoveringAnnotations(self):
        '''
        input:
            workspace = continuous workspace of size 11100 * 100
            annotations = covering the complete workspace, of size 100, 1000, and 10000 and alternating
            segments = segments of size 100, randomly distributed in workspace
        '''

        sizes = (100, 1000, 10000)
        workspace_size = sum(sizes) * 100
        segment_size = 1
        segment_spacing = 1000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        # workspace of size 1000000
        workspaces.add("default", "chr1",
                       SegmentList(iter=[(0, workspace_size)],
                                   normalize=True))

        # segments every 100bp
        segments.add("default", "chr1",
                     SegmentList(iter=[(x, x + segment_size) for x in range(0, workspace_size, segment_spacing)],
                                 normalize=True))

        x = 0
        intervals = collections.defaultdict(list)
        while x < workspace_size:
            for i, y in enumerate(sizes):
                intervals[i].append((x, x + y))
                x += y

        for i, y in enumerate(sizes):
            annotations.add("size-%i" % y, "chr1",
                            SegmentList(iter=intervals[i],
                                        normalize=True))

        self.check(workspaces["default"], annotations, segments)


class TestStatsTheAnnotator(TestStatsGat):

    def getSamples(self,
                   segments,
                   annotations,
                   workspace):

        fsegments, nsegments = tempfile.mkstemp()
        fworkspace, nworkspace = tempfile.mkstemp()
        fannotations, nannotations = tempfile.mkstemp()

        writeAnnotatorSegments2(fworkspace, workspace, "workspace")
        writeAnnotatorSegments2(fsegments, segments["default"], "segments")
        writeAnnotatorAnnotations(fannotations, annotations)

        statement = " ".join((ANNOTATOR_CMD,
                              "-seed %i" % random.getrandbits(16),
                              "-dumpsegments",
                              "-iterations %i" % self.sample_size,
                              "-annotation %s" % nannotations,
                              "-segments %s" % nsegments,
                              # increment to avoid removing results with P=1.0
                              "-baseline 1.01",
                              "-workspace %s" % nworkspace))

        process = subprocess.Popen(statement,
                                   cwd=os.getcwd(),
                                   shell=True,
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # process.stdin.close()
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise OSError("Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" %
                          (-process.returncode, stderr, statement))

        annotator_results = []

        intersection = SegmentList(clone=segments["default"]["chr1"])
        intersection.intersect(workspace["chr1"])
        segment_size = intersection.sum()

        # intersection = SegmentList( clone = annotations["default"]["chr1"] )
        # intersection.intersect( workspace["chr1"] )
        # annotation_size = intersection.sum()
        # print "segment_size", segment_size, annotation_size

        keep = False
        for line in stdout.split("\n"):

            if line.startswith("#"):
                continue
            if line.startswith("Z"):
                keep = True
                continue
            if not keep:
                continue

            data = line.split("\t")
            if len(data) != 9:
                break
            z, fold, pvalue, observed, expected, lower95, upper95, stddev, annotation = data
            pvalue, observed, expected, lower95, upper95, stddev = \
                list(map(float, (pvalue, observed, expected, lower95, upper95, stddev)))

            try:
                fold = observed / expected
            except ZeroDivisionError:
                fold = 0.0

            r = gat.DummyAnnotatorResult()
            r.track = "default"
            r.observed = observed * segment_size
            r.expected = expected * segment_size
            r.annotation = annotation[1:]
            r.fold = fold
            r.lower95 = lower95 * segment_size
            r.upper95 = upper95 * segment_size
            r.pvalue = pvalue
            r.qvalue = 1.0
            r.stddev = stddev * segment_size

            annotator_results.append(r)

        os.unlink(nsegments)
        os.unlink(nworkspace)
        os.unlink(nannotations)

        return annotator_results


class TestTrimming(GatTest):

    def testEndTrim(self):
        '''test end trimming.

        5 segments of 10 with a gap of 10.

        Count trimmed residues
        '''

        nrepeats = 100
        npoints = 100
        segment = 10
        return

        for size in range(2, 3 * segment):

            plt.figure()
            full = sum([list(range(x + size, x + segment - size))
                        for x in range(0, npoints, 2 * segment)], [])
            trimmed = sum([list(range(x, min(npoints - 1, x + size))) for x in range(0, npoints, 2 * segment)] +
                          [list(range(max(0, x + segment - size), x + segment)) for x in range(0, npoints, 2 * segment)], [])

            total = numpy.zeros(len(trimmed), numpy.int)

            for repeat in range(nrepeats):
                counts = numpy.zeros(npoints, numpy.int)
                for y in range(0, npoints, 2 * segment):
                    for point in range(y, y + segment):
                        ss = [(x, x + segment)
                              for x in range(0, npoints, 2 * segment)]
                        s = SegmentList(iter=ss,
                                        normalize=True)
                        s.trim_ends(point, size, numpy.random.randint(0, 2))
                        for start, end in s:
                            counts[start:end] += 1
                t = counts[trimmed]
                f = counts[full]
                if len(f) > 0:
                    self.assertEqual(min(f), max(
                        f), "min=%i, npoints=%i, size=%i" % (min(f), npoints, size))
                total += t

            plt.plot(range(len(total)), total)

        plt.show()


class TestEnrichmentGat(GatTest):
    '''
    check if enrichment values are correct by using workspace covering annotations.
    '''

    # complexity is sample_size * ntests
    # number of samples per test
    sample_size = 10
    # number of tests
    ntests = 100

    counter = None

    def getSamples(self,
                   segments,
                   annotations,
                   workspace):

        workspace_size = workspace["chr1"].sum()

        sampler = GatEngine.SamplerAnnotator(
            bucket_size=1, nbuckets=workspace_size)

        if self.counter == "NucleotideOverlap":
            counter = GatEngine.CounterNucleotideOverlap()
        elif self.counter == "SegmentOverlap":
            counter = GatEngine.CounterSegmentOverlap()
        counters = [counter]

        result = []

        workspace_generator = GatEngine.UnconditionalWorkspace()

        for x in range(self.ntests):
            result.append(
                gat.run(segments,
                        annotations,
                        workspace,
                        sampler,
                        counters,
                        workspace_generator,
                        num_samples=self.sample_size))

        return result

    def check(self, workspace, annotations, segments):

        if self.counter == None:
            return

        workspace_size = workspace["chr1"].sum()

        # annotator_results is two-level list
        # first level: samples
        # second level: annotations
        annotator_results = self.getSamples(segments,
                                            annotations,
                                            workspace)

        # make into a single list
        annotator_results = GatEngine.IOTools.flatten(annotator_results)

        tmp = SegmentList(clone=segments["default"]["chr1"])
        tmp.intersect(workspace["chr1"])
        segment_size = tmp.sum()

        annos = sorted(annotations.keys())
        expected, observed = {}, {}

        annotation_sizes = [annotations[x]["chr1"].sum() for x in annos]
        min_size = float(min(annotation_sizes))
        total_size = sum(annotation_sizes)
        nsegments = len(segments["default"]["chr1"])

        # collect observed and expected / scale by annotation size
        for annotation in annos:

            annotation_size = float(annotations[annotation]["chr1"].sum())
            nannotations = len(annotations[annotation]["chr1"])

            r = [x for x in annotator_results if x.annotation == annotation]
            assert len(r) > 0

            # proportion of total workspace
            expected_nucleotide_overlap = segment_size * annotation_size / total_size

            # bernoulli: E() = np
            expected_segment_overlap = min(nannotations, float(
                nannotations) * float(annotation_size) / total_size)
            avg = numpy.mean([x.expected for x in r])

            if self.counter == "NucleotideOverlap":
                scale = expected_nucleotide_overlap
            elif self.counter == "SegmentOverlap":
                scale = expected_segment_overlap

            expected[annotation] = [x.expected / scale for x in r]
            observed[annotation] = r[0].observed / scale

        plt.figure()
        plt.boxplot([expected[x] for x in annos])
        spacing = 1.0 / len(annos)
        width = spacing * 0.2
        for y, annotation in enumerate(annos):
            yy = spacing / 2.0 + y * spacing
            plt.axhline(observed[annotation], yy -
                        width, yy + width, color='g')
            plt.axhline(1.0, 0, 1, color='r', ls='--')

        all_vals = list(observed.values())
        for x in list(expected.values()):
            all_vals.extend(x)
        # plt.ylim( min(all_vals) - 10, max(all_vals) + 10 )
        if self.counter == "NucleotideOverlap":
            plt.yscale('log')

        filename = getPlotFilename(str(self))
        plt.savefig(filename)

        exp = sum([numpy.mean(expected[x]) for x in annos])
        obs = sum(observed[x] for x in annos)
        d = abs(obs - exp) / float(exp)
        self.assertTrue(d < 0.1, "variation too large :  %f / %f = %f > 10%%" %
                     (obs,
                      exp,
                      d))

    def checkUniformSegments(self,
                             annotation_sizes,
                             annotation_gap,
                             workspace_size,
                             segment_size,
                             segment_spacing):
        '''helper function for creating uniformly distributed segments.

        annotation_sizes - size of annotations
        workspace - fully covers annotations up to workspace_size
        '''

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        x = 0
        workspace = SegmentList()
        intervals = collections.defaultdict(list)
        while x < workspace_size:
            for i, y in enumerate(annotation_sizes):
                intervals[i].append((x, x + y))
                x += y + annotation_gap
                workspace.add(x, x + y)

        workspace.normalize()
        workspaces.add("default", "chr1", workspace)

        segments.add("default", "chr1",
                     SegmentList(iter=[(x, x + segment_size) for x in range(0, workspace_size, segment_spacing)],
                                 normalize=True))

        for i, y in enumerate(annotation_sizes):
            annotations.add("anno-%i" % i, "chr1",
                            SegmentList(iter=intervals[i],
                                        normalize=True))

        self.check(workspaces["default"],
                   annotations,
                   segments)

    def testVariableSizedAnnotationsWithoutEnrichment(self):
        '''
            workspace = annotations * 100
            annotations = 100, 1000, 10000, 20000, 30000
            segments = segments of size 100, every 1000 bases
        '''

        sizes = (100, 200, 500, 1000, 2000, 5000,
                 10000, 20000, 50000)  # = 500000)
        annotation_gap = 0
        workspace_size = 100 * sum(sizes)
        segment_size = 100
        segment_spacing = 1000

        self.checkUniformSegments(
            sizes, annotation_gap, workspace_size, segment_size, segment_spacing)

    def testEqualSizedAnnotationsWithoutEnrichment(self):
        '''
            workspace = continuous workspace of size 100000
            annotations = 4 * size of 1000
            segments = segments of size 100, every 1000 bases
        '''

        sizes = (1000, 1000, 1000, 1000)
        workspace_size = 100 * sum(sizes)
        annotation_gap = 0
        segment_size = 1
        segment_spacing = 1000

        self.checkUniformSegments(
            sizes, annotation_gap, workspace_size, segment_size, segment_spacing)

    def testEqualSizedAnnotationsWithGappedWorkspace(self):

        sizes = (1000, 1000, 1000, 1000)
        annotation_gap = 1000
        segment_size = 1
        segment_spacing = 1000
        workspace_size = 100 * sum(sizes)

        self.checkUniformSegments(
            sizes, annotation_gap, workspace_size, segment_size, segment_spacing)

    def testEqualSizedAnnotationsWithEnrichment(self):
        '''
            workspace = continuous workspace of size 100000
            annotations = 4 * size of 1000
            segments = size 100, enrichment in annotations 0 and 2
        '''

        sizes = (1000, 1000, 1000, 1000)
        annotation_gap = 0
        workspace_size = sum(sizes) * 100
        segment_size = 100
        segment_spacing = 1000

        workspaces, segments, annotations = \
            GatEngine.IntervalCollection( "workspace" ), \
            GatEngine.IntervalCollection( "segment" ), \
            GatEngine.IntervalCollection("annotation")

        workspaces.add("default", "chr1",
                       SegmentList(iter=[(0, workspace_size)],
                                   normalize=True))

        # add segments in first half of workspace as before
        intervals = [(x, x + segment_size)
                     for x in range(0, workspace_size / 2, segment_spacing)]
        # only add overlap with segments 0 and 2
        intervals.extend([(x, x + segment_size)
                          for x in range(workspace_size / 2, workspace_size, segment_spacing * 2)])
        segments.add("default", "chr1",
                     SegmentList(iter=intervals,
                                 normalize=True))

        x = 0
        intervals = collections.defaultdict(list)
        while x < workspace_size:
            for i, y in enumerate(sizes):
                intervals[i].append((x, x + y))
                x += y

        for i, y in enumerate(sizes):
            annotations.add("anno-%i" % i, "chr1",
                            SegmentList(iter=intervals[i],
                                        normalize=True))

        self.check(workspaces["default"],
                   annotations,
                   segments)


class TestEnrichmentGatNucleotideOverlap(TestEnrichmentGat):
    counter = "NucleotideOverlap"


class TestEnrichmentGatSegmentOverlap(TestEnrichmentGat):
    counter = "SegmentOverlap"

TestData = collections.namedtuple(
    "testdata", "segment_size, segment_spacing, annotation_size1, annotation_size2, annotation_gap_between, annotation_gap_within, annotation_mult")


class TestSpacingGat(GatTest):
    '''test relative sizes of segments versus annotations.

    If the annotations cover the workspace completely, the enrichment in overlap
    in one annotation must be offset by the depletion in another annotation.

    This test iterates over a large range of possible values varying:

    1. the segment sizes
    2. the size of the annotations
    3. the arrangement of the annotations (single, versus several consecutive ones)
    '''

    counter = "NucleotideOverlap"

    # number of samples calling sampler
    nsamples_sampler = 10

    # number of times pvalue is checked
    nsamples_pvalues = 100

    # number of samples of sampler
    nsamples = 10

    # total size of workspace
    workspace_size = 10000

    def build(self, p):
        '''build regular alternating annotations.'''

        intervals = [list() for x in range(2)]
        nsizes = 2
        sizes = [p.annotation_size1, p.annotation_size2]
        workspace = SegmentList()
        start = 0

        while start < self.workspace_size:
            for c, d in enumerate(sizes):
                for i in range(p.annotation_mult):
                    intervals[c].append((start, start + d))
                    workspace.add(start, start + d)
                    start += d + p.annotation_gap_within
                start -= p.annotation_gap_within
                start += p.annotation_gap_between

        # merges adjacent segments in workspace
        workspace.merge()

        # get effective workspace size
        workspace_size = workspace.max()

        segments = SegmentList(iter=[(x, x + p.segment_size)
                                     for x in range(0, workspace_size, p.segment_spacing)],
                               normalize=True)

        annotations = {}

        for c, y in enumerate(sizes):
            annotations["anno-%i" % c] = SegmentList(iter=intervals[c],
                                                     normalize=True)
        return segments, workspace, annotations

    def getSamples(self,
                   segments,
                   workspace,
                   annotations,
                   nsamples):

        sampler = GatEngine.SamplerAnnotator(bucket_size=1, nbuckets=100000)

        if self.counter == "NucleotideOverlap":
            counter = GatEngine.CounterNucleotideOverlap()
        elif self.counter == "SegmentOverlap":
            counter = GatEngine.CounterSegmentOverlap()

        counters = [counter]
        workspace_generator = GatEngine.UnconditionalWorkspace()

        result = collections.defaultdict(list)
        for x in range(self.nsamples_sampler):
            r = gat.run(segments,
                        annotations,
                        workspace,
                        sampler,
                        counters,
                        workspace_generator=workspace_generator,
                        num_samples=nsamples)

            for anno, v in r["default"].items():
                result[anno].append(v)

        return result

    def runSampler(self, segments, workspace, annotations, nsamples):

        ss = GatEngine.IntervalCollection("segment")
        ss.add("default", "chr1", segments)
        aa = GatEngine.IntervalCollection("annotation")
        for key, x in annotations.items():
            aa.add(key, "chr1", x)

        ww = GatEngine.IntervalCollection("workspace")
        ww.add("default", "chr1", workspace)

        samples = self.getSamples(ss, ww["default"], aa, nsamples)

        return samples

    def checkCoverage(self, params, outf):

        segments, workspace, annotations = self.build(params)

        ntries = 0

        nsamples = self.nsamples

        expected_ratio = float(
            annotations["anno-0"].sum()) / annotations["anno-1"].sum()

        while ntries < 2:

            samples = self.runSampler(
                segments, workspace, annotations, nsamples)

            observed, expected = [], []
            for anno in sorted(samples.keys()):
                rr = samples[anno]
                expected.append(numpy.mean([x.expected for x in rr]))
                observed.append(numpy.mean([x.observed for x in rr]))

            total_observed = sum(observed)
            total_expected = sum(expected)

            r1 = getDiff(total_expected, total_observed)

            observed_ratio = float(expected[0]) / expected[1]
            r2 = getDiff(expected_ratio, observed_ratio)

            if r2 > 0.1:
                if ntries == 1:
                    print("segs", segments)
                    print("work", workspace)
                    print("anno-0", annotations["anno-0"])
                    print("anno-1", annotations["anno-1"])
                    sys.exit(0)

                ntries += 1
                nsamples *= 10
            else:
                outf.write("\t".join(map(str,  ("\t".join(map(str, params)),
                                                total_observed,
                                                total_expected,
                                                r1, r2,
                                                observed, expected))) + "\n")
                outf.flush()
                break

    def checkPValues(self, params, outf):

        segments, workspace, annotations = self.build(params)

        # get randomly placed segments
        random_segments = self.runSampler(segments, workspace, annotations,
                                          self.nsamples_pvalues)

        print(random_segments["anno-1"][0])

    def testAnnotations(self):
        # segment_size, segment_spacing,
        segment_spacing = 10

        # params = TestData._make( (1,
        #                           10,
        #                           1,
        #                           2,
        #                           0,
        #                           0,
        #                           1,))
        # self.check( params )

        outf = open(getFilename(str(self)) + ".tsv", "w")
        outf.write("\t".join(TestData._fields + ("observed", "expected", "total_diff", "anno_diff",
                                                 "observed_vals", "expected_vals")) + "\n")

        for segment_size in (1, 10, 100):
            min_size = max(1, segment_size // 2)
            for annotation_size1 in (min_size, segment_size, segment_size * 2):
                for annotation_size2 in (annotation_size1, annotation_size1 * 2):
                    for annotation_gap_between in (0, min_size, segment_size, segment_size * 2):
                        for annotation_gap_within in (0, segment_size // 2, segment_size, segment_size * 2):

                            for annotation_mult in (1, 2, 3):
                                params = TestData._make((segment_size,
                                                         segment_spacing,
                                                         annotation_size1,
                                                         annotation_size2,
                                                         annotation_gap_within,
                                                         annotation_gap_between,
                                                         annotation_mult,))
                                # self.checkCoverage( params, outf )
                                self.checkPValues(params, outf)

        outf.close()

# class TestSpacingTheAnnotator( TestSpacingGat ):

#     def getSamples( self,
#                     segments,
#                     workspace,
#                     annotations,
#                     nsamples ):

#         result = collections.defaultdict( list )
#         for x in xrange( self.nsamples_sampler ):
#             r = getAnnotatorSamples( segments, annotations, workspace, nsamples )

#             for anno, v in r["default"].iteritems():
#                 result[anno].append( v )

#         return result


if __name__ == '__main__':
    unittest.main()
