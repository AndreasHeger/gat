import re
import os
import glob
import gat
import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats as Stats
import numpy

import gat.SegmentList as SegmentList
import gat.Engine as Engine

try:
    import matplotlib.pyplot as plt
    HASPLOT = True
except (ImportError, RuntimeError):
    HASPLOT = False


def dumpStats(coll, section, options):
    if section in options.output_stats or \
            "all" in options.output_stats or \
            len([x for x in options.output_stats
                 if re.search(x, section)]) > 0:
        coll.outputStats(E.openOutputFile(section))


def dumpBed(coll, section, options):
    if section in options.output_bed or \
            "all" in options.output_bed or \
            len([x for x in options.output_bed if re.search(x, section)]) > 0:
        coll.save(E.openOutputFile(section + ".bed"))


def readSegmentList(label,
                    filenames,
                    enable_split_tracks=False,
                    ignore_tracks=False):
    """read one or more segment files.

    Arguments
    ---------
    label : string
        Label to use for IntervalCollection.
    filenames : list
        List of filenames to load in :term:`bed` format.
    enable_split_tracks : bool
        If True, allow tracks to be split across multiple files.
    ignore_tracks : int
        If True, ignore track information.

    Returns
    -------
    segments : IntervalCollection
        The segment collection.
    """
    results = Engine.IntervalCollection(name=label)
    E.info("%s: reading tracks from %i files" % (label, len(filenames)))
    results.load(filenames,
                 allow_multiple=enable_split_tracks,
                 ignore_tracks=ignore_tracks)
    E.info("%s: read %i tracks from %i files" %
           (label, len(results), len(filenames)))
    return results


def readAnnotatorResults(filename):
    '''load annotator results from a tab-separated results table.'''

    annotator_results = []

    with IOTools.openFile(filename, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            r = gat.DummyAnnotatorResult._fromLine(line)
            annotator_results.append(r)

    return annotator_results


def expandGlobs(infiles):
    return IOTools.flatten([glob.glob(x) for x in infiles])


def buildSegments(options):
    '''load segments, annotations and workspace from parameters
    defined in *options*.

    The workspace will be split by isochores.

    returns segments, annotations and workspace.
    '''

    options.segment_files = expandGlobs(options.segment_files)
    options.annotation_files = expandGlobs(options.annotation_files)
    options.workspace_files = expandGlobs(options.workspace_files)
    options.sample_files = expandGlobs(options.sample_files)

    ##################################################
    # arguments sanity check
    if not options.segment_files:
        raise ValueError("please specify at least one segment file")
    if not options.annotation_files:
        raise ValueError("please specify at least one annotation file")
    if not options.workspace_files:
        raise ValueError("please specify at least one workspace file")

    # read one or more segment files
    segments = readSegmentList("segments",
                               options.segment_files,
                               ignore_tracks=options.ignore_segment_tracks)
    segments.normalize()

    if segments.sum() == 0:
        E.critical("no segments in input file - run aborted")
        raise ValueError("segments file is empty - run aborted")

    if len(segments) > 1000:
        raise ValueError(
            "too many (%i) segment files - use track definitions "
            "or --ignore-segment-tracks" % len(segments))

    annotations = readSegmentList(
        "annotations", options.annotation_files,
        enable_split_tracks=options.enable_split_tracks,
        ignore_tracks=options.annotations_label is not None)

    if options.annotations_label is not None:
        annotations.setName(options.annotations_label)

    if options.annotations_to_points:
        annotations.toPositions(options.annotations_to_points)

    if options.overlapping_annotations:
        # only sort, do not merge
        annotations.sort()
    else:
        annotations.normalize()

    workspaces = readSegmentList(
        "workspaces", options.workspace_files, options,
        options.enable_split_tracks)
    workspaces.normalize()

    # intersect workspaces to build a single workspace
    E.info("collapsing workspaces")
    dumpStats(workspaces, "stats_workspaces_input", options)
    workspaces.collapse()
    dumpStats(workspaces, "stats_workspaces_collapsed", options)

    # use merged workspace only, discard others
    workspaces.restrict("collapsed")

    # build isochores or intersect annotations/segments with workspace
    if options.isochore_files:

        # read one or more isochore files
        isochores = Engine.IntervalCollection(name="isochores")
        E.info("%s: reading isochores from %i files" %
               ("isochores", len(options.isochore_files)))
        isochores.load(options.isochore_files)
        dumpStats(isochores, "stats_isochores_raw", options)

        # merge isochores and check if consistent (fully normalized)
        isochores.sort()

        # check that there are no overlapping segments within isochores
        isochores.check()

        # TODO: flag is_normalized not properly set
        isochores.normalize()

        # check that there are no overlapping segments between isochores

        # truncate isochores to workspace
        # crucial if isochores are larger than workspace.
        isochores.intersect(workspaces["collapsed"])

    else:
        isochores = None

    return segments, annotations, workspaces, isochores


def applyIsochores(segments, annotations, workspaces,
                   options,
                   isochores=None,
                   truncate_segments_to_workspace=False,
                   truncate_workspace_to_annotations=False,
                   restrict_workspace=False,
                   ):
    '''apply isochores to segments and annotations.

    Segments and annotations are filtered in place to keep only those
    overlapping the workspace.

    If *isochores* are given, isochores are applied.

    If *truncate_segments_to_workspace*, truncate segments
    to workspace.

    If *restrict_workspace* is set, the workspace is confined
    to those parts that overlap both a segment and an annotation.

    If *truncate_workspace_to_annotations* is set, the workspace
    is truncated to keep only those parts that overlap annotations.

    returns a workspace divided into isochores.

    '''

    if isochores:
        # intersect isochores and workspaces, segments and annotations
        # workspace and annotations are truncated
        # with segments it is optional.
        E.info("adding isochores to workspace")
        workspaces.toIsochores(isochores, truncate=True)
        annotations.toIsochores(isochores, truncate=True)
        segments.toIsochores(
            isochores, truncate=options.truncate_segments_to_workspace)

        if workspaces.sum() == 0:
            raise ValueError("isochores and workspaces do not overlap")
        if annotations.sum() == 0:
            raise ValueError("isochores and annotations do not overlap")
        if segments.sum() == 0:
            raise ValueError("isochores and segments do not overlap")

        dumpStats(workspaces, "stats_workspaces_isochores", options)
        dumpStats(annotations, "stats_annotations_isochores", options)
        dumpStats(segments, "stats_segments_isochores", options)

        dumpBed(workspaces, "workspaces_isochores", options)
        dumpBed(annotations, "annotations_isochores", options)
        dumpBed(segments, "segments_isochores", options)

    else:
        # intersect workspace and segments/annotations
        # annotations and segments are truncated by workspace
        if options.truncate_segments_to_workspace:
            segments.intersect(workspaces["collapsed"])
        else:
            segments.filter(workspaces["collapsed"])

        annotations.intersect(workspaces["collapsed"])

        dumpStats(annotations, "stats_annotations_truncated", options)
        dumpStats(segments, "stats_segments_truncated", options)

    workspace = workspaces["collapsed"]

    if restrict_workspace:

        E.info("restricting workspace")
        # this is very cumbersome - refactor merge and collapse
        # to return an IntervalDictionary instead of adding it
        # to the list of tracks
        for x in (segments, annotations):
            if "merged" in segments:
                workspace.filter(segments["merged"])
            else:
                segments.merge()
                workspace.filter(segments["merged"])
                del segments["merged"]

        dumpStats(workspaces, "stats_workspaces_restricted", options)

    if truncate_workspace_to_annotations:

        E.info("truncating workspace to annotations")
        annotations.merge()
        annotations["merged"].normalize()
        workspace.intersect(annotations["merged"])
        del annotations["merged"]

        dumpStats(workspaces, "stats_workspaces_truncated", options)

    # segments.dump( open("segments_dump.bed", "w" ) )
    # workspaces.dump( open("workspaces_dump.bed", "w" ) )

    # output overlap stats
    # output segment densities per workspace
    if "overlap" in options.output_stats or \
            "all" in options.output_stats:
        for track in segments.tracks:
            workspaces.outputOverlapStats(
                E.openOutputFile("overlap_%s" % track),
                segments[track])

    return workspace


def readDescriptions(options):
    '''read descriptions from tab separated file.'''

    description_header, descriptions, description_width = [], {}, 0
    if options.input_filename_descriptions:
        E.info("reading descriptions from %s" %
               options.input_filename_descriptions)

        with IOTools.openFile(options.input_filename_descriptions) as inf:
            first = True
            for line in inf:
                if line.startswith("#"):
                    continue
                data = line[:-1].split("\t")

                if description_width:
                    assert len(data) - 1 == description_width, \
                        "inconsistent number of descriptions in %s" %\
                        options.input_filename_descriptions
                else:
                    description_width = len(data) - 1

                if first:
                    description_header = data[1:]
                    first = False
                else:
                    descriptions[data[0]] = data[1:]
        assert len(description_header) == description_width, \
            "number of descriptions (%i) inconsistent with header (%s) in %s" % \
            (description_width, len(description_header),
             options.input_filename_descriptions)

    return description_header, descriptions, description_width


class SegmentsSummary:

    '''summarize segments in a workspace.
    '''
    header = ("all_segments",
              "all_nucleotides",
              "segments_overlapping_workspace",
              "nucleotides_overlapping_workspace",
              "segments_outside_workspace",
              "nucleotides_outside_workspace",
              "truncated_segments",
              "truncated_nucleotides",
              "density_workspace",
              "proportion_truncated_segments",
              "proportion_extending_nucleotides",
              "summary_all_segments",
              "summary_segments_overlapping_workspace",
              "summary_truncated_segments")

    def __init__(self):
        pass

    def update(self, segments, workspace):
        '''compute summary statistics between two segments lists *segments* and *workspace*.
        This method computes:
            * number of segments/nucleotides
            * number of segments/nucleotides overlapping workspace
            * number of segments split by a boundary
            * number of nucleotides extending outside boundary
            * segment size distribution untruncated at boundaries (min, max, mean, media, q1, q3)
            * segment size distributino truncated at boundaries
        '''

        self.all_segments = len(segments)
        self.all_nucleotides = segments.sum()

        # build segments overlapping workspace
        segments_overlapping_workspace = SegmentList.SegmentList(
            clone=segments)
        segments_overlapping_workspace.filter(workspace)

        # build segments truncated by workspace
        truncated_segments = SegmentList.SegmentList(
            clone=segments_overlapping_workspace)
        truncated_segments.intersect(workspace)

        segments_extending_workspace = SegmentList.SegmentList(
            clone=segments)
        segments_extending_workspace.subtract(truncated_segments)

        # compute numbers
        self.segments_overlapping_workspace = len(truncated_segments)
        self.nucleotides_overlapping_workspace = truncated_segments.sum()

        self.segments_outside_workspace = self.all_segments - \
            self.segments_overlapping_workspace
        self.nucleotides_outside_workspace = self.all_nucleotides - \
            self.nucleotides_overlapping_workspace

        self.truncated_segments = len(segments_extending_workspace)
        self.truncated_nucleotides = segments_extending_workspace.sum()

        self.summary_all_segments = segments.summarize()
        self.summary_segments_overlapping_workspace = segments_overlapping_workspace.summarize()
        self.summary_truncated_segments = truncated_segments.summarize()

        workspace_size = workspace.sum()

        self.density_workspace, self.proportion_truncated_segments, self.proportion_extending_nucleotides = 0, 0, 0
        # some quality metrics
        if workspace_size > 0:
            self.density_workspace = float(
                self.nucleotides_overlapping_workspace) / workspace_size
        if self.segments_overlapping_workspace > 0:
            self.proportion_truncated_segments = float(
                self.truncated_segments) / self.segments_overlapping_workspace
            self.proportion_extending_nucleotides = float(
                self.truncated_nucleotides) / segments_overlapping_workspace.sum()

    def __str__(self):
        return "\t".join(["%i" % self.all_segments,
                          "%i" % self.all_nucleotides,
                          "%i" % self.segments_overlapping_workspace,
                          "%i" % self.nucleotides_overlapping_workspace,
                          "%i" % self.segments_outside_workspace,
                          "%i" % self.nucleotides_outside_workspace,
                          "%i" % self.truncated_segments,
                          "%i" % self.truncated_nucleotides,
                          "%f" % self.density_workspace,
                          "%f" % self.proportion_truncated_segments,
                          "%f" % self.proportion_extending_nucleotides,
                          "%s" % str(self.summary_all_segments),
                          "%s" % str(
                              self.summary_segments_overlapping_workspace),
                          "%s" % str(self.summary_truncated_segments)])


def outputMetrics(outfile, segments, workspace, track, section):
    '''output summary metrics

    Outputs summary metrics for segments in workspace. 
    .'''

    stats_per_isochore = []
    for isochore, ss in segments.items():
        stats = SegmentsSummary()
        stats.update(ss, workspace[isochore])
        stats_per_isochore.append(stats)

    for attribute in ("all_segments", "all_nucleotides",
                      "segments_overlapping_workspace",
                      "nucleotides_overlapping_workspace",
                      "nucleotides_outside_workspace",
                      "truncated_segments",
                      "truncated_nucleotides",
                      "density_workspace",
                      "proportion_truncated_segments",
                      "proportion_extending_nucleotides",
                      ):
        values = [getattr(x, attribute) for x in stats_per_isochore]
        outfile.write("%s\t%s\t%s\t%s\n" %
                      (track, section, attribute, Stats.Summary(values)))

    outfile.flush()


def outputResults(results,
                  options,
                  header,
                  description_header,
                  description_width,
                  descriptions,
                  format_observed="%i"):
    '''compute FDR and output results.'''

    pvalues = [x.pvalue for x in results]

    ##################################################
    ##################################################
    ##################################################
    # compute global fdr
    ##################################################
    E.info("computing FDR statistics")
    qvalues = Engine.getQValues(pvalues,
                                   method=options.qvalue_method,
                                   vlambda=options.qvalue_lambda,
                                   pi0_method=options.qvalue_pi0_method)

    try:
        results = [x._replace(qvalue=qvalue)
                   for x, qvalue in zip(results, qvalues)]
        is_tuple = True
    except AttributeError:
        # not a namedtuple
        for x, qvalue in zip(results, qvalues):
            x.qvalue = qvalue
            x.format_observed = format_observed

        is_tuple = False

    counters = set([x.counter for x in results])

    for counter in counters:

        if len(counters) == 1:
            outfile = options.stdout
            output = results
        else:
            outfilename = re.sub("%s", counter, options.output_tables_pattern)
            E.info("output for counter %s goes to outfile %s" %
                   (counter, outfilename))
            outfile = IOTools.openFile(outfilename, "w")
            output = [x for x in results if x.counter == counter]

        outfile.write(
            "\t".join(list(header) + list(description_header)) + "\n")

        if options.output_order == "track":
            output.sort(key=lambda x: (x.track, x.annotation))
        elif options.output_order == "observed":
            output.sort(key=lambda x: x.observed)
        elif options.output_order == "annotation":
            output.sort(key=lambda x: (x.annotation, x.track))
        elif options.output_order == "fold":
            output.sort(key=lambda x: x.fold)
        elif options.output_order == "pvalue":
            output.sort(key=lambda x: x.pvalue)
        elif options.output_order == "qvalue":
            output.sort(key=lambda x: x.qvalue)
        else:
            raise ValueError("unknown sort order %s" % options.output_order)

        for result in output:
            if is_tuple:
                outfile.write("\t".join(map(str, result)))
            else:
                outfile.write(str(result))

            if descriptions:
                try:
                    outfile.write(
                        "\t" + "\t".join(descriptions[result.annotation]))
                except KeyError:
                    outfile.write("\t" + "\t".join([""] * description_width))
            outfile.write("\n")

        if outfile != options.stdout:
            outfile.close()


def plotResults(results, options):
    '''plot annotator results.'''

    ##################################################
    # plot histograms
    if options.output_plots_pattern and HASPLOT:

        def buildPlotFilename(options, key):
            filename = re.sub("%s", key, options.output_plots_pattern)
            filename = re.sub("[^a-zA-Z0-9-_./]", "_", filename)
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)
            return filename

        E.info("plotting sample stats")

        for r in results:

            plt.figure()
            k = []
            if r.track != "merged":
                k.append(r.track)
            k.append(r.annotation)
            if r.counter != "na":
                k.append(r.counter)
            key = "-".join(k)

            s = r.samples
            hist, bins = numpy.histogram(s,
                                         bins=100)

            # plot bars
            plt.hist(s, bins=100, normed=True, label=key)

            plt.axvline(r.observed, color='r', linewidth=2)

            # plot estimated
            sigma = r.stddev
            mu = r.expected
            plt.plot(bins,
                     1.0 / (sigma * numpy.sqrt(2 * numpy.pi)) *
                     numpy.exp(- (bins - mu) ** 2 / (2 * sigma ** 2)),
                     label="std distribution",
                     linewidth=2,
                     color='g')

            plt.legend()
            filename = buildPlotFilename(options, key)
            plt.savefig(filename)

        E.info("plotting P-value distribution")

        key = "pvalue"
        plt.figure()

        x, bins, y = plt.hist([r.pvalue for r in results],
                              bins=numpy.arange(0, 1.05, 0.025),
                              label="pvalue")

        plt.hist([r.qvalue for r in results],
                 bins=numpy.arange(0, 1.05, 0.025),
                 label="qvalue",
                 alpha=0.5)

        plt.legend()

        # hist, bins = numpy.histogram( \
        #     [r.pvalue for r in Engine.iterator_results(annotator_results) ],
        #     bins = 20 )
        # plt.plot( bins[:-1], hist, label = key )

        filename = buildPlotFilename(options, key)
        plt.savefig(filename)
