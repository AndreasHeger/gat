#!/usr/bin/env python
##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
great - run great analysis
==========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares one or more genomic segments of interest 
against one more other genomic annotations.

GREAT says that a segment overlaps an annotation if the midpoint
of the segment overlaps the annotation.

The GREAT analysis is done on a per-isochore basis and overall.
For the overall analysis, there is no isochore correction
as the union of all isochores is treated as a single workspace.

This script implies several counters (``--counter`` option):

binom
    binomial model (as implemented by GREAT). Returns the probability
    of a certain number of segments overlapping (mid-point overlap) 
    an annotation

hyperg
    hypergeometric model. Returns the probabibility of a certain
    number of nucleotides within segments overlapping an annotation.
    Conceptually, all genomic positions are treated independently and 
    segments are fragmented into a collection of single bases that are
    placed indepently. This is valid for point intervals such as SNPs,
    but of course ignores segment size, though it might be used as a
    guide towards the expected fold enrichment, that can be computed
    more accurately using :mod:`gat-run.py`

Usage
-----

Example::

   python gat-great.py 
      --segment-file=segments.bed.gz 
      --workspace-file=workspace.bed.gz 
      --annotation-file=annotations_architecture.bed.gz  
 
Type::

   python gat-run.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import collections
import types
import glob
import time
import numpy
import numpy.random

import gat
import gat.Experiment as E
import gat.IOTools as IOTools
import gat.IO as IO
import GatSegmentList

import scipy.stats

GREAT_RESULT = collections.namedtuple("GREAT",
                                      ("track",
                                       "annotation",
                                       "isochore",
                                       "counter",
                                       # nsegments_overlapping_annotations
                                       "observed",
                                       "expected",
                                       "nsegments_in_workspace",
                                       "nannotations_in_workspace",
                                       "nsegments_overlapping_annotation",
                                       "nannotations_overlapping_segments",
                                       "basecoverage_intersection",
                                       "basecoverage_segments",
                                       "basecoverage_annotation",
                                       "basecoverage_workspace",
                                       "fraction_coverage_annotation",
                                       "fold",
                                       "pvalue",
                                       "qvalue"))


def main(argv):

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-a", "--annotation-file", "--annotations", dest="annotation_files", type="string", action="append",
                      help="filename with annotations [default=%default].")

    parser.add_option("-s", "--segment-file", "--segments", dest="segment_files", type="string", action="append",
                      help="filename with segments. Also accepts a glob in parentheses [default=%default].")

    parser.add_option("-w", "--workspace-file", "--workspace", dest="workspace_files", type="string", action="append",
                      help="filename with workspace segments. Also accepts a glob in parentheses [default=%default].")

    parser.add_option("-i", "--isochore-file", "--isochores", dest="isochore_files", type="string", action="append",
                      help="filename with isochore segments. Also accepts a glob in parentheses [default=%default].")

    parser.add_option("-o", "--order", dest="output_order", type="choice",
                      choices=(
                          "track", "annotation", "fold", "pvalue", "qvalue"),
                      help="order results in output by fold, track, etc. [default=%default].")

    parser.add_option("-q", "--qvalue-method", dest="qvalue_method", type="choice",
                      choices=(
                          "storey", "BH", "bonferroni", "holm", "hommel", "hochberg", "BY", "none"),
                      help="method to perform multiple testing correction by controlling the fdr [default=%default].")

    parser.add_option("--qvalue-lambda", dest="qvalue_lambda", type="float",
                      help="fdr computation: lambda [default=%default].")

    parser.add_option("--qvalue-pi0-method", dest="qvalue_pi0_method", type="choice",
                      choices=("smoother", "bootstrap"),
                      help="fdr computation: method for estimating pi0 [default=%default].")
    parser.add_option("--descriptions", dest="input_filename_descriptions", type="string",
                      help="filename mapping annotation terms to descriptions. "
                      " if given, the output table will contain additional columns "
                      " [default=%default]")

    parser.add_option("--ignore-segment-tracks", dest="ignore_segment_tracks", action="store_true",
                      help="ignore segment tracks - all segments belong to one track [default=%default]")

    parser.add_option("--enable-split-tracks", dest="enable_split_tracks", action="store_true",
                      help="permit the same track to be in multiple files [default=%default]")

    parser.add_option("--output-bed", dest="output_bed", type="choice", action="append",
                      choices=("all",
                               "annotations", "segments",
                               "workspaces", "isochores",
                               "overlap"),
                      help="output bed files [default=%default].")

    parser.add_option("--output-stats", dest="output_stats", type="choice", action="append",
                      choices=("all",
                               "annotations", "segments",
                               "workspaces", "isochores",
                               "overlap"),
                      help="output overlap summary stats [default=%default].")

    parser.add_option("--restrict-workspace", dest="restrict_workspace", action="store_true",
                      help="restrict workspace to those segments that contain both track"
                      " and annotations [default=%default]")

    parser.add_option("--counter", dest="counters", type="choice", action="append",
                      choices=("binom", "hyperg"),
                      help="counter to use [default=%default].")

    parser.add_option("--output-tables-pattern", dest="output_tables_pattern", type="string",
                      help="output pattern for result tables. Used if there are multiple counters used [default=%default].")

    parser.set_defaults(
        annotation_files=[],
        segment_files=[],
        workspace_files=[],
        sample_files=[],
        counters=[],
        output_stats=[],
        output_bed=[],
        output_tables_pattern="%s.tsv.gz",
        output_order="fold",
        input_filename_counts=None,
        input_filename_results=None,
        pvalue_method="empirical",
        output_plots_pattern=None,
        output_samples_pattern=None,
        qvalue_method="storey",
        qvalue_lambda=None,
        qvalue_pi0_method="smoother",
        ignore_segment_tracks=False,
        input_filename_descriptions=None,
        conditional="unconditional",
        conditional_extension=None,
        conditional_expansion=None,
        restrict_workspace=False,
        enable_split_tracks=False,
        shift_expansion=2.0,
        shift_extension=0,
        overlap_mode="midpoint",
        truncate_workspace_to_annotations=False,
        truncate_segments_to_workspace=False
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    tstart = time.time()

    if len(options.counters) == 0:
        options.counters.append("binom")

    ############################################
    segments, annotations, workspaces, isochores = IO.buildSegments(options)
    E.info("intervals loaded in %i seconds" % (time.time() - tstart))

    # filter segments by workspace
    workspace = IO.applyIsochores(
        segments, annotations, workspaces, options, isochores)

    ############################################
    description_header, descriptions, description_width = IO.readDescriptions(
        options)

    ############################################
    ############################################
    # compute per contig

    # compute bases covered by workspace
    workspace2basecoverage, isochores = {}, []
    for contig, ww in workspace.iteritems():
        workspace2basecoverage[contig] = ww.sum()
        isochores.append(contig)

    # compute percentage of bases covered by annotations in workspace
    # per isochore
    annotation2basecoverage = collections.defaultdict(dict)
    for annotation, aa in annotations.iteritems():
        for isochore, a in aa.iteritems():
            # need to truncate to workspace?
            annotation2basecoverage[annotation][isochore] = a.sum()

    results_per_contig = collections.defaultdict(list)

    E.info("computing counts per isochore")
    # results per isochore

    def emptyResult(segment, annotation, isochore,
                    counter,
                    nsegments_in_workspace,
                    basecoverage_annotation,
                    basecoverage_workspace):
        return GREAT_RESULT._make((
            segment, annotation, isochore,
            counter,
            0,  # observed
            0,  # expected
            nsegments_in_workspace,
            0,  # nannotations_in_workspace
            0,  # nsegments_overlapping_annotation
            0,  # nannotations_overlapping_segments
            0,  # basecoverage_intersection
            0,  # basecoverage_segments
            basecoverage_annotation,
            basecoverage_workspace,
            0.0,
            1.0,
            1.0,
            1.0))

    for isochore in isochores:
        basecoverage_workspace = workspace2basecoverage[isochore]

        # iterate over all isochores
        for segment, segmentdict in segments.iteritems():
            try:
                ss = segmentdict[isochore]
                # select segments overlapping workspace
                segments_in_workspace = GatSegmentList.SegmentList(clone=ss)
                segments_in_workspace.intersect(workspace[isochore])
                # number of segments in workspace
                nsegments_in_workspace = len(segments_in_workspace)
            except KeyError:
                ss = None

            basecoverage_segments = segments_in_workspace.sum()

            for annotation, annotationdict in annotations.iteritems():

                # if annotation != "GO:0030957": continue

                try:
                    aa = annotationdict[isochore]
                except KeyError:
                    aa = None

                # p_A: proportion of bases covered by annotation
                try:
                    basecoverage_annotation = annotation2basecoverage[
                        annotation][isochore]
                except KeyError:
                    basecoverage_annotation = 0

                if ss == None or aa == None:
                    for counter in options.counters:
                        results_per_contig[(counter, segment, annotation)].append(emptyResult(segment, annotation,
                                                                                              isochore,
                                                                                              counter,
                                                                                              nsegments_in_workspace,
                                                                                              basecoverage_annotation,
                                                                                              basecoverage_workspace))
                    continue

                # select segments overlapping annotation
                segments_overlapping_annotation = GatSegmentList.SegmentList(
                    clone=ss)
                segments_overlapping_annotation.intersect(
                    annotations[annotation][isochore])
                # number of segments in annotation
                nsegments_overlapping_annotation = ss.intersectionWithSegments(annotations[annotation][isochore],
                                                                               mode=options.overlap_mode)

                # number of nucleotides at the intersection of segments,
                # annotation and workspace
                basecoverage_intersection = segments_overlapping_annotation.sum()

                annotations_overlapping_segments = GatSegmentList.SegmentList(
                    clone=aa)
                annotations_overlapping_segments.intersect(ss)
                nannotations_overlapping_segments = len(
                    annotations_overlapping_segments)

                nannotations_in_workspace = len(aa)
                if nannotations_in_workspace == 0:
                    for counter in options.counters:
                        results_per_contig[(counter, segment, annotation)].append(emptyResult(segment,
                                                                                              annotation,
                                                                                              isochore,
                                                                                              counter,
                                                                                              nsegments_in_workspace,
                                                                                              basecoverage_annotation,
                                                                                              basecoverage_workspace))
                    continue

                fraction_coverage_annotation = basecoverage_annotation / \
                    float(basecoverage_workspace)
                fraction_hit_annotation = float(
                    nannotations_overlapping_segments) / nannotations_in_workspace

                for counter in options.counters:
                    if counter.startswith("binom"):
                        # GREAT binomial probability over "regions"
                        # n = number of genomic regions = nannotations_in_workspace
                        # ppi = fraction of genome annotated by annotation = fraction_coverage_annotation
                        # kpi = genomic regions with annotation hit by segments = nannotations_in_segments
                        # sf = survival functions = 1 -cdf
                        # probability of observing >kpi in a sample of n where the probabily of succes is
                        # ppi.
                        pvalue = scipy.stats.binom.sf(nsegments_overlapping_annotation - 1,
                                                      nsegments_in_workspace,
                                                      fraction_coverage_annotation)

                        expected = fraction_coverage_annotation * \
                            nsegments_in_workspace
                        observed = nsegments_overlapping_annotation

                    elif counter.startswith("hyperg"):

                        # hypergeometric probability over nucleotides
                        # Sampling without replacement
                        # x,M,n,M
                        # x = observed number of nucleotides in overlap of segments,annotations and workspace
                        # M = number of nucleotides in workspace
                        # n = number of nucleotides in annotations (and workspace)
                        # N = number of nucleotides in segments (and workspace)
                        # P-value of obtaining >x number of nucleotides
                        # overlapping.
                        rv = scipy.stats.hypergeom(basecoverage_workspace,
                                                   basecoverage_annotation,
                                                   basecoverage_segments)

                        pvalue = rv.sf(basecoverage_intersection)
                        expected = rv.mean()
                        observed = basecoverage_intersection

                    if expected != 0:
                        fold = float(observed) / expected
                    else:
                        fold = 1.0

                    r = GREAT_RESULT._make((
                        segment, annotation, isochore,
                        counter,
                        observed,
                        expected,
                        nsegments_in_workspace,
                        nannotations_in_workspace,
                        nsegments_overlapping_annotation,
                        nannotations_overlapping_segments,
                        basecoverage_intersection,
                        basecoverage_segments,
                        basecoverage_annotation,
                        basecoverage_workspace,
                        fraction_coverage_annotation,
                        fold,
                        pvalue,
                        1.0))
                    # print "\t".join( map(str, r))
                    results_per_contig[
                        (counter, segment, annotation)].append(r)

    E.info("merging counts per isochore")

    # compute sums
    results = []

    for niteration, pair in enumerate(results_per_contig.iteritems()):

        counter, segment, annotation = pair[0]
        data = pair[1]

        nsegments_in_workspace = sum([x.nsegments_in_workspace for x in data])
        nsegments_overlapping_annotation = sum([x.observed for x in data])
        nannotations_in_workspace = sum(
            [x.nannotations_in_workspace for x in data])
        nannotations_overlapping_segments = sum(
            [x.nannotations_overlapping_segments for x in data])

        basecoverage_intersection = sum(
            [x.basecoverage_intersection for x in data])
        basecoverage_segments = sum([x.basecoverage_segments for x in data])
        basecoverage_annotation = sum(
            [x.basecoverage_annotation for x in data])
        basecoverage_workspace = sum([x.basecoverage_workspace for x in data])

        fraction_coverage_annotation = basecoverage_annotation / \
            float(basecoverage_workspace)

        if counter.startswith("binom"):
            pvalue = scipy.stats.binom.sf(nsegments_overlapping_annotation - 1,
                                          nsegments_in_workspace,
                                          fraction_coverage_annotation)
            expected = fraction_coverage_annotation * nsegments_in_workspace
            observed = nsegments_overlapping_annotation
        elif counter.startswith("hyperg"):
            rv = scipy.stats.hypergeom(basecoverage_workspace,
                                       basecoverage_annotation,
                                       basecoverage_segments)

            pvalue = rv.sf(basecoverage_intersection)
            expected = rv.mean()
            observed = basecoverage_intersection

        if expected != 0:
            fold = float(observed) / expected
        else:
            fold = 1.0

        r = GREAT_RESULT._make((
            segment, annotation, "all",
            counter,
            observed,
            expected,
            nsegments_in_workspace,
            nannotations_in_workspace,
            nsegments_overlapping_annotation,
            nannotations_overlapping_segments,
            basecoverage_intersection,
            basecoverage_segments,
            basecoverage_annotation,
            basecoverage_workspace,
            fraction_coverage_annotation,
            fold,
            pvalue,
            1.0))

        results.append(r)

    IO.outputResults(results,
                     options,
                     GREAT_RESULT._fields,
                     description_header,
                     description_width,
                     descriptions)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
