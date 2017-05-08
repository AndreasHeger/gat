import os
import re
import optparse
import collections
import gzip
import numpy

import gat.Bed as Bed
import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats as Stats
import gat.IO as IO

import gat.Engine as Engine

import multiprocessing.pool
import multiprocessing


def readFromBedOld(filenames, name="track"):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *contig* and *track*.

    If no track is given, the *name* attribute is taken.
    '''

    segment_lists = collections.defaultdict(
        lambda: collections.defaultdict(Engine.SegmentList))

    if name == "track":
        f = lambda x: x.mTrack["name"]
    elif name == "name":
        f = lambda x: x.mFields[0]
    else:
        raise ValueError("unknown name: '%s'" % name)

    for filename in filenames:
        infile = IOTools.openFile(filename, "r")
        for bed in Bed.iterator(infile):
            try:
                name = f(bed)
            except TypeError:
                name = "default"
            segment_lists[name][bed.contig].add(bed.start, bed.end)

    return segment_lists


class OptionGroup(optparse.OptionGroup):
    pass


def buildParser(usage=None):
    '''return gat command line parser.
    '''

    parser = optparse.OptionParser(version="%prog version: $Id:",
                                   usage=usage)

    group = OptionGroup(parser, "Input options")

    group.add_option(
        "-a", "--annotation-bed-file", "--annotations", "--annotation-file",
        dest="annotation_files", type="string", action="append",
        help="filename with annotations [default=%default].")

    group.add_option(
        "-s", "--segment-bed-file", "--segments", "--segment-file",
        dest="segment_files", type="string",
        action="append",
        help="filename with segments. Also accepts a "
        "glob in parentheses [default=%default].")

    group.add_option(
        "-w", "--workspace-bed-file", "--workspace", "--workspace-file",
        dest="workspace_files", type="string", action="append",
        help="filename with workspace segments. Also "
        "accepts a glob in parentheses [default=%default].")

    group.add_option(
        "-i", "--isochore-bed-file", "--isochores", "--isochore-file",
        dest="isochore_files", type="string", action="append",
        help="filename with isochore segments. Also "
        "accepts a glob in parentheses [default=%default].")

    group.add_option(
        "-l", "--sample-file", dest="sample_files",
        type="string", action="append",
        help="filename with sample files. Start processing "
        "from samples [default=%default].")

    group.add_option(
        "--input-counts-file", dest="input_filename_counts",
        type="string",
        help="start processing from counts - no segments "
        "required [default=%default].")

    group.add_option(
        "--input-results-file", dest="input_filename_results",
        type="string",
        help="start processing from results - no segments "
        "required [default=%default].")

    group.add_option(
        "--ignore-segment-tracks", dest="ignore_segment_tracks",
        action="store_true",
        help="ignore segment tracks - all segments belong "
        "to one track and called 'merged' [default=%default]")

    group.add_option(
        "--with-segment-tracks", dest="ignore_segment_tracks",
        action="store_false",
        help="the segments data file is arranged in tracks. "
        "[default=%default]")

    group.add_option(
        "--enable-split-tracks", dest="enable_split_tracks",
        action="store_true",
        help="permit the same track to be in multiple "
        "files [default=%default]")

    group.add_option(
        "--overlapping-annotations", dest="overlapping_annotations",
        action="store_true",
        help="the annotations within a track are overlapping and should not "
        "be merged. This is useful for working with short-read data. "
        "[default=default]")

    group.add_option(
        "--annotations-label", dest="annotations_label",
        type="string",
        help="ignore tracks in annotations and instead set them "
        "to label "
        "[default=default]")

    group.add_option(
        "--annotations-to-points", dest="annotations_to_points",
        type="choice",
        choices=("midpoint", "start", "end"),
        help="convert annotations from segments to positions. Available "
        "methods are 'midpoint', 'start' or 'end'. "
        "[default=default]")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Output options")

    group.add_option(
        "-o", "--order", dest="output_order", type="choice",
        choices=(
            "track", "annotation", "fold", "pvalue", "qvalue"),
        help="order results in output by fold, track, etc. "
        "[default=%default].")

    group.add_option(
        "--output-tables-pattern", dest="output_tables_pattern",
        type="string",
        help="output pattern for result tables. Used if there "
        "are multiple counters used [default=%default].")

    group.add_option(
        "--output-counts-pattern", dest="output_counts_pattern",
        type="string",
        help="output pattern for counts [default=%default].")

    group.add_option(
        "--output-plots-pattern", dest="output_plots_pattern",
        type="string",
        help="output pattern for plots [default=%default]")

    group.add_option(
        "--output-samples-pattern",
        dest="output_samples_pattern", type="string",
        help="output pattern for samples. Samples are "
        "stored in bed format, one for "
        " each segment [default=%default]")

    group.add_option(
        "--output-stats", dest="output_stats", type="choice",
        action="append",
        choices=("all",
                 "annotations", "segments",
                 "workspaces", "isochores",
                 "overlap",
                 "sample",
                 "segment_metrics",
                 "sample_metrics"),
        help="output overlap summary stats [default=%default].")

    group.add_option(
        "--output-bed", dest="output_bed", type="choice",
        action="append",
        choices=("all",
                 "annotations", "segments",
                 "workspaces", "isochores",
                 "overlap"),
        help="output bed files [default=%default].")

    group.add_option(
        "--descriptions", dest="input_filename_descriptions",
        type="string",
        help="filename mapping annotation terms to "
        "descriptions. "
        "if given, the output table will contain additional "
        "columns "
        "[default=%default]")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Sampling algorithm options")

    group.add_option(
        "-c", "--counter", dest="counters", type="choice",
        action="append",
        choices=("nucleotide-overlap",
                 "nucleotide-density",
                 "segment-overlap",
                 "segment-midoverlap",
                 "annotation-overlap",
                 "annotation-midoverlap"),
        help="quantity to use for estimating enrichment "
        "[default=%default].")

    group.add_option(
        "-m", "--sampler", dest="sampler", type="choice",
        choices=("annotator",
                 "segments",
                 "shift",
                 "local-permutation",
                 "global-permutation",
                 "uniform",
                 "brute-force"),
        help="quantity to test [default=%default].")

    group.add_option(
        "-n", "--num-samples", dest="num_samples", type="int",
        help="number of samples to compute [default=%default].")

    group.add_option(
        "--shift-extension", dest="shift_extension",
        type="float",
        help="if the sampling method is 'shift', create a "
        "segment of size # anound the segment "
        "to determine the size of the region for "
        "shifthing [default=%default].")

    group.add_option(
        "--shift-expansion", dest="shift_expansion",
        type="float",
        help="if the sampling method is 'shift', multiply each "
        "segment by # "
        "to determine the size of the region for "
        "shifthing [default=%default].")

    group.add_option(
        "--bucket-size", dest="bucket_size", type="int",
        help="size of a bin for histogram of segment lengths. "
        "If 0, it will be automatically "
        "scaled to fit nbuckets [default=%default]")

    group.add_option(
        "--nbuckets", dest="nbuckets", type="int",
        help="number of bins for histogram of segment "
        "lengths [default=%default]")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Statistics options")

    group.add_option(
        "-p", "--pvalue-method", dest="pvalue_method",
        type="choice",
        choices=("empirical", "norm", ),
        help="type of pvalue reported [default=%default].")

    group.add_option(
        "-q", "--qvalue-method", dest="qvalue_method",
        type="choice",
        choices=(
            "storey", "BH", "bonferroni", "holm", "hommel",
            "hochberg", "BY", "none"),
        help="method to perform multiple testing correction "
        "by controlling the fdr [default=%default].")

    group.add_option(
        "--qvalue-lambda", dest="qvalue_lambda", type="float",
        help="fdr computation: lambda [default=%default].")

    group.add_option(
        "--qvalue-pi0-method", dest="qvalue_pi0_method",
        type="choice",
        choices=("smoother", "bootstrap"),
        help="fdr computation: method for estimating pi0 "
        "[default=%default].")

    group.add_option(
        "--pseudo-count", dest="pseudo_count", type="float",
        help="pseudo count. The pseudo count is added to both "
        "the observed and expected overlap. "
        "Using a pseudo-count avoids gat reporting fold changes "
        "of 0 [default=%default].")

    group.add_option(
        "--null", dest="null", type="string",
        help="null hypothesis. The default is to test "
        "categories "
        "for enrichment/depletion. "
        "If a filename with gat output is given, gat will test "
        "for the difference "
        "in fold change between the segments supplied and in "
        "the other file [default=%default].")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Processing options")

    group.add_option(
        "-e", "--cache", dest="cache", type="string",
        help="filename for caching samples [default=%default].")

    group.add_option(
        "-t", "--num-threads", dest="num_threads", type="int",
        help="number of threads to use for sampling "
        "[default=%default]")

    group.add_option(
        "--random-seed", dest='random_seed', type="int",
        help="random seed to initialize number generator "
        "with [%default].")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Workspace manipulation (experimental)")

    group.add_option(
        "--conditional", dest="conditional", type="choice",
        choices=("unconditional", "annotation-centered",
                 "segment-centered", "cooccurance"),
        help="conditional workspace creation [default=%default]"
        "*cooccurance* - compute enrichment only within "
        "workspace "
        "segments that contain both segments "
        "and annotations, "
        "*annotation-centered* - workspace centered around "
        "annotations. See --conditional-extension,"
        "segment-centered - workspace centered around "
        "segments. See --conditional-extension")

    group.add_option(
        "--conditional-extension", dest="conditional_extension",
        type="int",
        help="if workspace is created conditional, extend by "
        "this amount (in bp) [default=%default].")

    group.add_option(
        "--conditional-expansion", dest="conditional_expansion",
        type="float",
        help="if workspace is created conditional, expand by "
        "this amount (ratio) [default=%default].")

    group.add_option(
        "--restrict-workspace", dest="restrict_workspace",
        action="store_true",
        help="restrict workspace to those segments that "
        "contain both track "
        "and annotations [default=%default]")

    group.add_option(
        "--truncate-workspace-to-annotations",
        dest="truncate_workspace_to_annotations",
        action="store_true",
        help="truncate workspace with annotations "
        "[default=%default]")

    group.add_option(
        "--truncate-segments-to-workspace",
        dest="truncate_segments_to_workspace",
        action="store_true",
        help="truncate segments to workspace before "
        "sampling [default=%default]")

    parser.add_option_group(group)

    parser.set_defaults(
        annotation_files=[],
        annotations_label=None,
        annotations_to_points=None,
        bucket_size=0,
        cache=None,
        conditional="unconditional",
        conditional_expansion=None,
        conditional_extension=None,
        counters=[],
        enable_split_tracks=False,
        ignore_segment_tracks=True,
        input_filename_counts=None,
        input_filename_descriptions=None,
        input_filename_results=None,
        nbuckets=100000,
        null="default",
        num_samples=1000,
        num_threads=0,
        output_bed=[],
        output_counts_pattern=None,
        output_order="fold",
        output_plots_pattern=None,
        output_samples_pattern=None,
        output_stats=[],
        output_tables_pattern="%s.tsv.gz",
        overlapping_annotations=False,
        pseudo_count=1.0,
        pvalue_method="empirical",
        qvalue_lambda=None,
        qvalue_method="BH",
        qvalue_pi0_method="smoother",
        random_seed=None,
        restrict_workspace=False,
        sample_files=[],
        sampler="annotator",
        segment_files=[],
        shift_expansion=2.0,
        shift_extension=0,
        truncate_segments_to_workspace=False,
        truncate_workspace_to_annotations=False,
        workspace_files=[],
    )

    return parser


def iterator_results(annotator_results):
    '''iterate over all results.'''
    for k1, v1 in annotator_results.items():
        for k2, v2 in v1.items():
            yield v2


class DummyAnnotatorResult:

    format_observed = "%i"
    format_expected = "%6.4f"
    format_fold = "%6.4f"
    format_pvalue = "%6.4e"

    def __init__(self):
        pass

    @classmethod
    def _fromLine(cls, line):
        x = cls()
        data = line[:-1].split("\t")
        x.track, x.annotation = data[:2]
        x.counter = "na"
        x.observed, x.expected, x.lower95, x.upper95, x.stddev, x.fold, x.l2fold, x.pvalue, x.qvalue = \
            list(map(float, data[2:11]))
        if len(data) > 11:
            (track_nsegments,
             track_size,
             track_density,
             annotation_nsegments,
             annotation_size,
             annotation_density,
             overlap_nsegments,
             overlap_size,
             overlap_density,
             percent_overlap_nsegments_track,
             percent_overlap_size_track,
             percent_overlap_nsegments_annotation,
             percent_overlap_size_annotation) = list(map(float, data[11:24]))

        return x

    def __str__(self):
        return "\t".join((self.track,
                          self.annotation,
                          self.format_observed % self.observed,
                          self.format_expected % self.expected,
                          self.format_expected % self.lower95,
                          self.format_expected % self.upper95,
                          self.format_expected % self.stddev,
                          self.format_fold % self.fold,
                          self.format_pvalue % self.pvalue,
                          self.format_pvalue % self.qvalue))


WorkData = collections.namedtuple("WorkData",
                                  "track sample_id sampler segments "
                                  "annotations contig_annotations "
                                  "workspace contig_workspace "
                                  "counters ")


def computeSample(args):
    '''compute a single sample.
    '''

    workdata, samples_outfile, metrics_outfile, lock = args

    (track,
     sample_id,
     sampler,
     segs,
     annotations,
     contig_annotations,
     workspace,
     contig_workspace,
     counters) = workdata

    # E.debug("track=%s, sample=%s - started" % (track, str(sample_id)))

    counts = E.Counter()

    sample_id = str(sample_id)

    outf_samples = samples_outfile

    if samples_outfile:
        if lock:
            lock.acquire()
            outf_samples = IOTools.openFile(samples_outfile, "a")

        samples_outfile.write("track name=%s\n" % sample_id)

        if lock:
            outf_samples.close()
            lock.release()

    sample = Engine.IntervalDictionary()

    for isochore in list(segs.keys()):

        counts.pairs += 1

        # skip empty isochores
        if workspace[isochore].isEmpty or segs[isochore].isEmpty:
            counts.skipped += 1
            continue

        counts.sampled += 1
        r = sampler.sample(segs[isochore], workspace[isochore])

        # TODO : activate
        # self.outputSampleStats( sample_id, isochore, r )

        sample.add(isochore, r)

        # save sample
        if samples_outfile:
            if lock:
                lock.acquire()
                outf_samples = IOTools.openFile(samples_outfile, "a")

            for start, end in r:
                outf_samples.write("%s\t%i\t%i\n" % (isochore, start, end))

            if lock:
                outf_samples.close()
                lock.release()

    # re-combine isochores
    # adjacent intervals are merged.
    sample.fromIsochores()

    if metrics_outfile:
        if lock:
            lock.acquire()
            outf = IOTools.openFile(metrics_outfile, "a")
        else:
            outf = metrics_outfile

        IO.outputMetrics(outf, sample, workspace, track, sample_id)

        if lock:
            outf.close()
            lock.release()

    counts_per_track = [collections.defaultdict(float) for x in counters]
    # compute counts for each counter
    for counter_id, counter in enumerate(counters):
        # TODO: choose aggregator
        for annotation in annotations.tracks:
            counts_per_track[counter_id][annotation] = sum([
                counter(sample[contig],
                        contig_annotations[annotation][contig],
                        contig_workspace[contig])
                for contig in list(sample.keys())])

    # E.debug("track=%s, sample=%s - completed" % (track,str(sample_id )))

    return counts_per_track


class UnconditionalSampler:

    def __init__(self,
                 num_samples,
                 samples,
                 samples_outfile,
                 sampler,
                 workspace_generator,
                 counters,
                 outfiles,
                 num_threads=1):
        self.num_samples = num_samples
        self.samples = samples
        self.samples_outfile = samples_outfile
        self.sampler = sampler
        self.workspace_generator = workspace_generator
        self.counters = counters

        self.outfile_sample_stats = outfiles.get("sample_stats", None)
        self.outfile_sample_metrics = outfiles.get("sample_metrics", None)

        if self.outfile_sample_stats:
            E.debug("sample stats go to %s" % self.outfile_sample_stats)
            self.outfile_sample_stats.write(
                "sample\tisochore\tnsegments\tnnucleotides\tmean\t"
                "std\tmin\tq1\tmedian\tq3\tmax\n")

        self.last_sample_id = None
        self.all_lengths = []
        self.num_threads = num_threads

    def outputSampleStats(self, sample_id, isochore, sample):

        def _write(sample_id, isochore, lengths):
            if self.outfile_sample_stats:
                self.outfile_sample_stats.write("\t".join(map(str, (
                    sample_id,
                    isochore,
                    len(lengths),
                    numpy.sum(lengths),
                    numpy.mean(lengths),
                    numpy.std(lengths),
                    min(lengths),
                    Stats.percentile(lengths, 0.25),
                    numpy.median(lengths),
                    Stats.percentile(lengths, 0.75),
                    max(lengths)))) + "\n")

        if sample_id != self.last_sample_id:
            if self.last_sample_id is not None:
                _write(self.last_sample_id, "all", numpy.sort(
                    numpy.array(self.all_lengths)))
            self.all_lengths = []
            self.last_sample_id = sample_id

        if sample_id:
            l = sample.asLengths()
            self.all_lengths.extend(l)
            _write(sample_id, isochore, numpy.sort(numpy.array(l)))

    def computeSamples(self, work, report_interval=100):
        '''compute samples according to work.

        returns a list of results.
        '''
        n = len(work)

        E.debug('sampling will work on %i items' % n)

        results = []

        if self.num_threads == 0:
            for i, w in enumerate(work):
                r = computeSample(
                    (w, self.samples_outfile, self.outfile_sample_metrics,
                     None))
                if i % report_interval == 0:
                    E.info("%i/%i done (%5.2f)" % (i, n, 100.0 * i / n))
                results.append(r)
        else:
            E.info("generating processpool with %i threads for %i items" %
                   (self.num_threads, len(work)))

            manager = multiprocessing.Manager()

            lock = manager.Lock()

            pool = multiprocessing.Pool(self.num_threads)

            # use file names - not files when multiprocessing
            samples_outfile, metrics_outfile = None, None
            if self.samples_outfile:
                samples_outfile = self.samples_outfile.name
                self.samples_outfile.flush()
            if self.outfile_sample_metrics:
                metrics_outfile = self.outfile_sample_metrics.name
                self.outfile_sample_metrics.flush()

            ww = [(w, samples_outfile, metrics_outfile, lock) for w in work]

            for i, r in enumerate(pool.imap_unordered(computeSample, ww)):
                if i % report_interval == 0:
                    E.info("%i/%i done (%5.2f)" % (i, n, 100.0 * i / n))
                results.append(r)

            pool.close()
            pool.join()

        return results

    def sample(self, track, counts, counters, segs,
               annotations, workspace,
               outfiles):
        '''sample and return counts.

        Return a list of counted results for each counter.
        '''

        E.info("performing unconditional sampling")
        counts_per_track = [collections.defaultdict(list) for x in counters]

        # rebuild non-isochore annotations and workspace
        contig_annotations = annotations.clone()
        contig_annotations.fromIsochores()
        contig_annotations.setName("contig_" + annotations.getName())

        contig_workspace = workspace.clone()
        contig_workspace.fromIsochores()

        E.info("workspace without conditioning: %i segments, %i nucleotides" %
               (workspace.counts(),
                workspace.sum()))

        temp_segs, _, temp_workspace = self.workspace_generator(
            segs, None, workspace)

        E.info("workspace after conditioning: %i segments, %i nucleotides" %
               (workspace.counts(),
                workspace.sum()))

        if workspace.sum() == 0:
            E.warn("empty workspace - no computation performed")
            return None

        work = [WorkData(track,
                         x,
                         self.sampler,
                         temp_segs,
                         annotations,
                         contig_annotations,
                         temp_workspace,
                         contig_workspace,
                         counters,
                         ) for x in range(self.num_samples)]

        if self.num_threads > 0:
            E.info("setting up shared data for multi-processing")
            annotations.share()
            contig_annotations.share()
            contig_workspace.share("contig_workspace")
            temp_segs.share("generated_segments")
            temp_workspace.share("generated_workspace")

        E.info("sampling started")
        results = self.computeSamples(work)
        E.info("sampling completed")

        if self.num_threads > 0:
            E.info("retrieving private data")
            annotations.unshare()
            contig_annotations.unshare()
            contig_workspace.unshare()
            temp_segs.unshare()
            temp_workspace.unshare()

        # collate results
        for result in results:
            for counter_id, counter in enumerate(counters):
                for annotation in annotations.tracks:
                    counts_per_track[counter_id][annotation].append(
                        result[counter_id][annotation])

        self.outputSampleStats(None, "", [])

        return counts_per_track


class ConditionalSampler(UnconditionalSampler):

    def sample(self, track, counts, counters, segs, annotations, workspace,
               outfiles):
        '''conditional sampling - sample using only those
        segments that contain both a segment and an annotation.

        return dictionary with counts per track
        '''

        E.info("performing conditional sampling")
        counts_per_track = [collections.defaultdict(list) for x in counters]

        # rebuild non-isochore annotations and workspace
        contig_annotations = annotations.clone()
        contig_annotations.fromIsochores()
        contig_annotations.setName("contig_" + annotations.getName())

        contig_workspace = workspace.clone()
        contig_workspace.fromIsochores()

        E.info("setting up shared data for multi-processing")
        annotations.share()
        contig_annotations.share()
        contig_workspace.share("contig_workspace")

        E.info("workspace without conditioning: %i segments, %i nucleotides" %
               (workspace.counts(),
                workspace.sum()))

        if workspace.sum() == 0:
            E.warn("empty workspace - no computation performed")
            return None

        # compute samples conditionally - need to proceed by annotation
        for annoid, annotation in enumerate(annotations.tracks):

            annos = annotations[annotation]

            temp_segs, temp_annotations, temp_workspace = \
                self.workspace_generator(segs, annos, workspace)

            # set up sharing
            temp_segs.share("generated_segments")
            temp_workspace.share("generated_workspace")

            E.info("workspace for annotation %s: %i segments, %i nucleotides" %
                   (annotation,
                    temp_workspace.counts(),
                    temp_workspace.sum()))

            work = [WorkData('_'.join((track, annoid)),
                             x,
                             self.sampler,
                             temp_segs,
                             annotations,
                             contig_annotations,
                             temp_workspace,
                             contig_workspace,
                             counters,
                             ) for x in range(self.num_samples)]

            E.info("sampling for annotation '%s' started" % annotation)
            results = self.computeSamples(work)
            E.info("sampling for annotation '%s' completed" % annotation)

            for result in results:
                for counter_id, counter in enumerate(counters):
                    counts_per_track[counter_id][annotation].append(
                        result[counter_id][annotation])

        return counts_per_track


def run(segments,
        annotations,
        workspace,
        sampler,
        counters,
        workspace_generator,
        **kwargs):
    '''run an enrichment analysis.

    segments: an IntervalCollection
    workspace: an IntervalCollection
    annotations: an IntervalCollection

    kwargs recognized are:

    cache
       filename of cache

    num_samples
       number of samples to compute

    output_counts_pattern
       output counts to filename

    output_samples_pattern
       if given, output samles to these files, one per segment

    sample_files
       if given, read samples from these files.

    fdr
       method to compute qvalues

    outfiles
       dictionary of optional additional output files.

    pseudo_count
       pseudo_count to add to observed and expected values

    reference
       data with reference observed and expected values.
    '''

    # get arguments
    num_samples = kwargs.get("num_samples", 10000)
    cache = kwargs.get("cache", None)
    output_counts_pattern = kwargs.get("output_counts_pattern", None)
    sample_files = kwargs.get("sample_files", [])
    pseudo_count = kwargs.get("pseudo_count", 1.0)
    reference = kwargs.get("reference", None)
    output_samples_pattern = kwargs.get("output_samples_pattern", None)
    outfiles = kwargs.get("outfiles", {})
    num_threads = kwargs.get("num_threads", 0)

    ##################################################
    ##################################################
    ##################################################
    # computing summary metrics for segments
    if "segment_metrics" in outfiles:
        E.info("computing summary metrics for segments")
        outfile = outfiles["segment_metrics"]
        outfile.write("track\tsection\tmetric\t%s\n" %
                      "\t".join(Stats.Summary().getHeaders()))
        for track in segments.tracks:
            IO.outputMetrics(outfile,
                             segments[track],
                             workspace,
                             track,
                             'segments',
                             )
        E.info("wrote summary metrics for segments to %s" % str(outfile))

    ##################################################
    ##################################################
    ##################################################
    # collect observed counts from segments
    E.info("collecting observed counts")
    observed_counts = []
    for counter in counters:
        observed_counts.append(Engine.computeCounts(
            counter=counter,
            aggregator=sum,
            segments=segments,
            annotations=annotations,
            workspace=workspace,
            workspace_generator=workspace_generator))

    ##################################################
    ##################################################
    ##################################################
    # sample and collect counts
    ##################################################
    E.info("starting sampling")

    if cache:
        E.info("samples are cached in %s" % cache)
        samples = Engine.SamplesCached(filename=cache)
    elif sample_files:
        if not output_samples_pattern:
            raise ValueError(
                "require output_samples_pattern if loading samples from files")
        # build regex
        regex = re.compile(re.sub("%s", "(\S+)", output_samples_pattern))
        E.info("loading samples from %i files" % len(sample_files))
        samples = Engine.SamplesFile(
            filenames=sample_files,
            regex=regex)
    else:
        samples = Engine.Samples()

    sampled_counts = {}

    counts = E.Counter()

    ntracks = len(segments.tracks)

    for ntrack, track in enumerate(segments.tracks):

        segs = segments[track]

        E.info("sampling: %s: %i/%i" % (track, ntrack + 1, ntracks))

        if output_samples_pattern and not sample_files:
            filename = re.sub("%s", track, output_samples_pattern)
            E.debug("saving samples to %s" % filename)
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)
            if filename.endswith(".gz"):
                samples_outfile = gzip.open(filename, "w")
            else:
                samples_outfile = open(filename, "w")
        else:
            samples_outfile = None

        if workspace_generator.is_conditional:
            outer_sampler = ConditionalSampler(num_samples,
                                               samples,
                                               samples_outfile,
                                               sampler,
                                               workspace_generator,
                                               counters,
                                               outfiles,
                                               num_threads=num_threads)
        else:
            outer_sampler = UnconditionalSampler(num_samples,
                                                 samples,
                                                 samples_outfile,
                                                 sampler,
                                                 workspace_generator,
                                                 counters,
                                                 outfiles,
                                                 num_threads=num_threads)

        counts_per_track = outer_sampler.sample(
            track, counts, counters, segs, annotations, workspace, outfiles)

        # skip empty tracks
        if counts_per_track is None:
            continue

        if samples_outfile:
            samples_outfile.close()

        sampled_counts[track] = counts_per_track

        # old code, refactor into loop to save samples
        if 0:
            E.info("sampling stats: %s" % str(counts))
            if track not in samples:
                E.warn("no samples for track %s" % track)
                continue

            # clean up samples
            del samples[track]

    E.info("sampling finished")

    # build annotator results
    E.info("computing PValue statistics")

    annotator_results = list()
    counter_id = 0
    for counter, observed_count in zip(counters, observed_counts):
        for track, r in observed_count.items():
            for annotation, observed in r.items():
                temp_segs, temp_annos, temp_workspace = workspace_generator(
                    segments[track],
                    annotations[annotation],
                    workspace)

                # ignore empty results
                if temp_workspace.sum() == 0:
                    continue

                # if reference is given, p-value will indicate difference
                # The test that track and annotation are present is done
                # elsewhere
                if reference:
                    ref = reference[track][annotation]
                else:
                    ref = None

                annotator_results.append(Engine.AnnotatorResultExtended(
                    track=track,
                    annotation=annotation,
                    counter=counter.name,
                    observed=observed,
                    samples=sampled_counts[track][counter_id][annotation],
                    track_segments=temp_segs,
                    annotation_segments=temp_annos,
                    workspace=temp_workspace,
                    reference=ref,
                    pseudo_count=pseudo_count))
        counter_id += 1

    # dump (large) table with counts
    if output_counts_pattern:
        for counter in counters:
            name = counter.name
            filename = re.sub("%s", name, output_counts_pattern)

            E.info("writing counts to %s" % filename)
            output = [x for x in annotator_results if x.counter == name]
            outfile = IOTools.openFile(filename, "w")
            outfile.write("track\tannotation\tobserved\tcounts\n")

            for o in output:
                outfile.write("%s\t%s\t%i\t%s\n" %
                              (o.track, o.annotation,
                               o.observed,
                               ",".join(["%i" % x for x in o.samples])))

    return annotator_results


def fromCounts(filename):
    '''build annotator results from a tab-separated table
    with counts.'''

    annotator_results = []

    with IOTools.openFile(filename, "r") as infile:

        E.info("loading data")

        header = infile.readline()
        if not header == "track\tannotation\tobserved\tcounts\n":
            raise ValueError("%s not a counts file: got %s" % (infile, header))

        for line in infile:
            track, annotation, observed, counts = line[:-1].split("\t")
            samples = numpy.array(
                list(map(float, counts.split(","))), dtype=numpy.float)
            observed = float(observed)
            annotator_results.append(Engine.AnnotatorResult(
                track=track,
                annotation=annotation,
                counter="na",
                observed=observed,
                samples=samples))

    return annotator_results
