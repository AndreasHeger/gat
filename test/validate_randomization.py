"""test gat implementation."""

import random
import tempfile
import shutil
import os
import re
import gzip
import sys
import itertools
import optparse
import gat
import numpy
import math
import multiprocessing
import SVGdraw

import Experiment as E


class Runner:
    '''run gat or other simulator.'''

    ntests = 100

    def run(self, segments, workspace):
        samples = []

        for x in range(self.ntests):
            sample = self.getSample(segments, workspace)
            samples.append(sample)

        return samples


class RunnerGat(Runner):

    def __init__(self):
        self.sampler = gat.SamplerAnnotator()

    def getSample(self, segments, workspace):
        '''return a single sample'''
        return self.sampler.sample(segments, workspace)


class SegmentGenerator:
    '''generate a test set.


    returns two normalized segment lists for segment and workspace

    nsegments: number of segments
    segment_length: segment length
    workspace_nregions: number of workspace regions
    workspace_length: size of workspace region
    workspace_gap: gap between workspace regions
    '''

    headers = ["class", "nsegments", "segment_length",
               "workspace_regions", "workspace_length", "workspace_gap",
               "workspace_size"]

    def __init__(self, nsegments, segment_length, workspace_nregions, workspace_length, workspace_gap):

        self.nsegments = nsegments
        self.segment_length = segment_length
        self.workspace_nregions = workspace_nregions
        self.workspace_length = workspace_length
        self.workspace_gap = workspace_gap

    def __str__(self):
        return "\t".join(map(str, ("gapped_workspace",
                                   self.nsegments,
                                   self.segment_length,
                                   self.workspace_nregions,
                                   self.workspace_length,
                                   self.workspace_gap,
                                   self.workspace_length * self.workspace_nregions)))

    def createWorkspace(self):
        '''create a workspace.'''

        result = []
        start, end = self.workspace_gap, self.workspace_gap
        for x in range(self.workspace_nregions):
            end = start + self.workspace_length
            result.append((start, end))
            start = end + self.workspace_gap

        return result

    def createSet(self):
        '''create a test set.

        Segments are constructed starting from the first residue in
        the workspace and then added with one residue gap in-between
        to prevent them from merging. 

        Segments are NOT randomly distributed within the workspace.

        Segments can partially overlap the workspace, but gaps in
        the workspace are accounted for properly.

        Returns a list of segments and workspaces.
        '''

        workspace = self.createWorkspace()

        assert len(workspace) > 0

        workspace_idx = 0
        start = max(0, workspace[workspace_idx][0] - (self.segment_length - 1))
        segments = []

        for x in range(self.nsegments):
            end = start + self.segment_length
            segments.append((start, end))

            # add gap between segments
            end += 1

            if end > workspace[workspace_idx][1]:
                # segment extending beyond current workspace segment.
                # advance workspace
                workspace_idx += 1
                # continue advancing workspace until end of workspace is larger than end of segment
                # accounts for segments straddling gaps in workspace or
                # multiple workspace segments.
                while workspace_idx < len(workspace) and workspace[workspace_idx][1] < end:
                    workspace_idx += 1

                if workspace_idx == len(workspace):
                    # can't find any more workspaces, done
                    break

                start = max(workspace[workspace_idx - 1][1] + 1,
                            workspace[workspace_idx][0] - (self.segment_length - 1))
            else:
                start = end

        if len(segments) != self.nsegments:
            # print segments
            # print workspace
            # E.warn( "not enough space in workspace for %i segments" % (self.nsegments) )
            pass

        _segments = gat.SegmentList(iter=segments, normalize=True)
        _workspace = gat.SegmentList(iter=workspace, normalize=True)

        return _segments, _workspace


class Validator(object):
    '''validate samples with respect to workspace and segments.'''

    # stringency level: percent of true value
    stringency_level = 0.1

    def __init__(self, segments, workspace, ntests):
        self.segments = segments
        self.workspace = workspace
        self.ntests = ntests


class ValidatorNumSamples(Validator):

    headers = ("samples_ok", "nsamples")

    def validate(self, samples):
        return "\t".join(("%i" % (len(samples) == self.ntests),
                          "%i" % len(samples)))


class ValidatorSegmentLength(Validator):
    '''check length distribution of segments.'''

    headers = ("length_ok",)

    def validate(self, samples):

        # check segment lengths
        l = [x[1] - x[0] for x in self.segments]
        values_input = min(l), max(l), numpy.mean(l), numpy.std(l)

        fail = False

        for i, sample in enumerate(samples):
            l = [x[1] - x[0] for x in sample]
            values_sample = min(l), max(l), numpy.mean(l), numpy.std(l)

            for val, inp, samp in zip(("min", "max", "mean", "std"),
                                      values_input,
                                      values_sample):
                d = abs(inp - samp) / float(inp)

                # segment length distribution fails
                if d >= self.stringency_level:
                    fail = True
                    E.warn("segment length distribution in sample %i: expected %s (%f) != observed %s (%f)" %
                           (i, val, inp, val, samp))

                    break

            if fail:
                break
        else:
            fail = False

        return "\t".join(("%i" % (not fail)), )


def computeSegmentDensityProfile(workspace,
                                 samples):
    '''compute sample counts within workspace.

    returns array of size workspace with counts
    '''

    l = workspace.max()
    counts = numpy.zeros(l, numpy.int)
    ntests = len(samples)

    segment_sizes = []
    starts = numpy.zeros(l + 1, numpy.int)
    ends = numpy.zeros(l + 1, numpy.int)

    for sample_id, s in enumerate(samples):
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

    counts_within_workspace = []
    for start, end in workspace:
        counts_within_workspace.extend(counts[start:end])

    return numpy.array(counts_within_workspace, dtype=numpy.int), segment_sizes, starts, ends


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


class SamplePlot:

    linewidth = 5
    linespace = 1

    def __init__(self):
        pass

    def plot(self,
             filename,
             segments,
             workspace,
             samples):
        nlines = 2 + len(samples)

        height = nlines * self.linewidth * self.linespace

        width = workspace.max()

        print(height, width)

        root = SVGdraw.drawing()
        canvas = SVGdraw.svg((0, 0, width, heigh), "100%", "100%")

        root.setSVG(canvas)

        root.toXml(filename)


def plotCounts(filename,
               counts_within_workspace,
               segment_sizes,
               starts, ends,
               workspace,
               expected_coverage=None,
               density=0):

    import matplotlib.pyplot as plt

    l = workspace.max()

    if l > 10:
        dx = 10
    else:
        dx = 1

    newy = smooth(counts_within_workspace, window_len=dx)

    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplots_adjust(right=0.7)
    # plt.axes( [0.1,0.1,0.51,0.5] )

    plt.subplot("311")
    plt.plot(range(len(counts_within_workspace)),
             counts_within_workspace, '.', label="coverage")

    plt.plot(range(len(counts_within_workspace)), newy, '-',
             label="smooth (%i)" % dx)

    plt.title("%s : density = %6.4f" % (filename, density))
    # : nworkspaces=%i, sample_length=%i" % ( filename, len(workspace), len(samples) ) )
    plt.xlabel("position")
    plt.ylabel("counts")
    if expected_coverage:
        d = expected_coverage * 0.1
        plt.plot(range(len(counts_within_workspace)), [expected_coverage] * len(counts_within_workspace),
                 '-', label="expected")
        plt.plot(range(len(counts_within_workspace)), [expected_coverage - d] * len(counts_within_workspace),
                 '--')
        plt.plot(range(len(counts_within_workspace)), [expected_coverage + d] * len(counts_within_workspace),
                 '--')
    plt.legend(loc=(1.03, 0.2))

    plt.subplot("312")
    segment_sizes.sort()
    segment_sizes = list(zip(*segment_sizes))
    plt.plot(segment_sizes[0], label="stddev")
    plt.plot(segment_sizes[1], label="min")
    plt.plot(segment_sizes[2], label="max")
    plt.plot(segment_sizes[3], label="mean")
    plt.plot(segment_sizes[4], label="median")
    plt.legend(loc=(1.03, 0.2))
    plt.xlabel("sample")
    plt.ylabel("segment size")

    plt.subplot("313")
    plt.plot(starts, label="starts")
    plt.plot(ends, label="ends")
    plt.legend(loc=(1.03, 0.2))
    plt.xlabel("position")
    plt.ylabel("counts")

    # plt.savefig( filename )
    plt.show()


class ValidatorSegmentDistribution(Validator):

    headers = ("nucleotide_ok",
               "average_coverage_ok",
               "uniform_coverage_ok",
               "overlap_ok",
               "average_dev",
               "uniform_dev",
               "segment_density",
               "average_overlap",
               "expected_overlap",
               "average_coverage",
               "expected_coverage")

    def computeExpectedCoverageOld(self, samples):

        # expected coverage
        # number_of_samples * segment_length / workspace_length
        # modifier:
        # can be simplified

        # number of segments
        nsegments = float(len(self.segments))

        # number of nucleotides in segments overlapping the workspace
        tmp = gat.SegmentList(clone=self.workspace)
        tmp.intersect(self.segments)
        segment_overlap = tmp.sum()

        # density of segment nucleotides in workspace
        segment_density = float(segment_overlap) / self.workspace.sum()

        segment_length = segment_size / nsegments
        expected = len(samples) * segment_overlap / float(self.workspace.sum())
        #   float(workspace.sum()) / (workspace.sum() + segments.sum()  * len(workspace) )

        return nsegments, segment_overlap, segment_density, segment_length, expected

    def computeExpectedCoverage(self, samples):

        # expected coverage
        # number_of_samples * segment_length / workspace_length
        # modifier:
        # can be simplified

        # number of segments
        nsegments = float(len(self.segments))

        # number of nucleotides in segments overlapping the workspace
        tmp = gat.SegmentList(clone=self.workspace)
        tmp.intersect(self.segments)
        expected_overlap = tmp.sum()

        # for computing expected overlap use complete segments
        # as there is no truncation
        segment_size = self.segments.sum()

        # density of segment nucleotides in workspace
        segment_density = float(expected_overlap) / self.workspace.sum()

        # average length of segments within workspace
        expected_segment_length = segment_size / nsegments

        # expected coverage of segments
        expected_coverage = len(samples) * \
            expected_overlap / float(self.workspace.sum())
        #   float(workspace.sum()) / (workspace.sum() + segments.sum()  * len(workspace) )

        return (nsegments, segment_size, segment_density,
                expected_overlap,
                expected_segment_length,
                expected_coverage)

    def validate(self, samples):

        # filename = getPlotFilename()

        # compute expected coverage
        (nsegments, segment_size, segment_density,
         expected_overlap,
         expected_segment_length,
         expected_coverage) = self.computeExpectedCoverage(samples)

        # compute actual coverage counts
        counts_within_workspace, segment_sizes, starts, ends = computeSegmentDensityProfile(self.workspace,
                                                                                            samples)

        # plotCounts( None,
        #             counts_within_workspace,
        #             segment_sizes,
        #             starts, ends,
        #             self.workspace,
        #             expected_coverage = expected_coverage,
        #             density = segment_density )

        ##################################
        # check if each sample has the correct number of nucleotides
        nucleotide_ok = True

        sums = [x.sum() for x in samples]
        overlaps = []
        for x, s in enumerate(samples):
            tmp = gat.SegmentList(clone=self.workspace, normalize=True)
            s.normalize()
            tmp.intersect(s)
            ovl = tmp.sum()
            overlaps.append(ovl)
            if ovl != expected_overlap:
                nucleotide_ok = False
                E.warn("incorrect number of nucleotides in sample %i, got %i, expected %i, %s" %
                       (x, ovl, expected_overlap, samples[x]))

        ##################################
        # check if average overlap is ok
        overlap_ok = True
        average_overlap = numpy.mean(overlaps)
        average_overlap_d = abs(
            average_overlap - expected_overlap) / float(expected_overlap)
        if average_overlap_d >= self.stringency_level:
            overlap_ok = False
            E.warn("expected_overlap counts (%f) != sampled counts (%f)" % (expected_overlap,
                                                                            average_overlap))

        ##################################
        # check average coverage
        average_ok = True
        average_coverage = counts_within_workspace.mean()
        average_d = abs(average_coverage - expected_coverage) / \
            float(expected_coverage)

        if average_d >= self.stringency_level:
            average_ok = False
            E.warn("expected_coverage counts (%f) != sampled counts (%f)" % (expected_coverage,
                                                                             counts_within_workspace.mean()))

        # check for uniform coverage
        uniform_ok = True
        stddev = numpy.std(counts_within_workspace)
        uniform_d = stddev / float(expected_coverage)

        if uniform_d >= self.stringency_level:
            uniform_ok = False
            E.warn("coverage variation too large : stddev (%f) / %f = %f > %f" %
                   (stddev,
                    expected_coverage,
                    uniform_d,
                    self.stringency_level))

        return "\t".join(("%i" % nucleotide_ok,
                          "%i" % average_ok,
                          "%i" % uniform_ok,
                          "%i" % overlap_ok,
                          "%6.4f" % average_d,
                          "%6.4f" % uniform_d,
                          "%6.4f" % segment_density,
                          "%6.4f" % average_overlap,
                          "%6.4f" % expected_overlap,
                          "%6.4f" % average_coverage,
                          "%6.4f" % expected_coverage))


def runSimulation(data):

    lock, outfile, segmentor, runner, validators = data

    segments, workspaces = segmentor.createSet()

    line = []
    line.append(str(segmentor))

    # print segments
    # print workspaces

    # gat is not multiprocessing safe - TODO
    # instantiate here
    samples = runner().run(segments, workspaces)

    p = SamplePlot()
    p.plot(segments, workspaces, samples)

    for validator in validators:
        v = validator(segments, workspaces, runner.ntests)
        line.append(v.validate(samples))

    lock.acquire()
    # using outfile fails - file closed after first job
    # finishes -> file closed error
    sys.stdout.write("\t".join(line) + "\n")
    lock.release()


class Test:
    '''run a test.'''

    def __init__(self, runner, test_generator, validators):
        self.runner = runner
        self.test_generator = test_generator
        self.validators = validators

    def run(self, outfile, processors=1):

        tasks = []

        manager = multiprocessing.Manager()
        lock = manager.Lock()

        for segmentor in self.test_generator:
            headers = segmentor.headers
            tasks.append((lock, outfile, segmentor,
                          self.runner, self.validators))

        for v in self.validators:
            headers.extend(v.headers)
        outfile.write("%s\n" % "\t".join(headers))

        E.info("created %i tasks for %i workers" % (len(tasks), processors))

        if processors > 1:
            pool = multiprocessing.Pool(processors)
            pool.map(runSimulation, tasks)
        else:
            for task in tasks:
                runSimulation(task)


def test_segmented_workspaces():
    '''5 workspaces.'''

    # number of workspaces
    nworkspaces = (5,)

    # workspace_length
    workspace_length = range(10, 100, 20)

    # workspace gaps
    workspace_gap = range(0, 100, 20)

    # segment lengths
    segment_length = range(1, 100, 20)

    # number of segments
    nsegments = range(1, 50, 10)

    for _nworkspaces, _workspace_length, _workspace_gap, _segment_length, _nsegments in \
            itertools.product(nworkspaces, workspace_length, workspace_gap, segment_length, nsegments):
        yield SegmentGenerator(segment_length=_segment_length,
                               nsegments=_nsegments,
                               workspace_nregions=_nworkspaces,
                               workspace_length=_workspace_length,
                               workspace_gap=_workspace_gap)


def small_test_segmented_workspaces():
    '''5 workspaces.'''

    # number of workspaces
    nworkspaces = (5,)

    # workspace_length
    workspace_length = (10,)

    # workspace gaps
    workspace_gap = (60,)

    # segment lengths
    segment_length = (21,)

    # number of segments
    nsegments = (1,)

    for _nworkspaces, _workspace_length, _workspace_gap, _segment_length, _nsegments in \
            itertools.product(nworkspaces, workspace_length, workspace_gap, segment_length, nsegments):
        yield SegmentGenerator(segment_length=_segment_length,
                               nsegments=_nsegments,
                               workspace_nregions=_nworkspaces,
                               workspace_length=_workspace_length,
                               workspace_gap=_workspace_gap)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id$",
                                   usage=globals()["__doc__"])

    parser.add_option("-p", "--proc", dest="processors", type="int",
                      help="use # processors [%default]")

    parser.set_defaults(
        processors=1)

    options, args = E.Start(parser, argv=argv)

    t1 = Test(RunnerGat,
              small_test_segmented_workspaces(),
              [ValidatorNumSamples,
               ValidatorSegmentDistribution])

    t1.run(options.stdout,
           processors=options.processors)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
