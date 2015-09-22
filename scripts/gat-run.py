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
'''gat-run - run the genomic annotation tool
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares one or more genomic segments of interest against
one more other genomic annotations.

Usage
-----

Example::

   python gat-run.py
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
import time
import numpy
import random

import gat
import gat.Experiment as E
import gat.IO as IO
import gat.Stats as Stats
import gat.SegmentList as SegmentList
import gat.Engine as Engine


def fromSegments(options, args):
    '''run analysis from segment files.

    This is the most common use case.
    '''

    tstart = time.time()

    # build segments
    segments, annotations, workspaces, isochores = IO.buildSegments(options)

    E.info("intervals loaded in %i seconds" % (time.time() - tstart))

    # open various additional output files
    outfiles = {}
    for section in ("sample",
                    "segment_metrics",
                    "sample_metrics",
                    ):
        if section in options.output_stats or \
            "all" in options.output_stats or \
                len([x for x in options.output_stats
                     if re.search(x, "section")]) > 0:
            outfiles[section] = E.openOutputFile(section)

    if 'sample_metrics' in outfiles:
        outfiles['sample_metrics'].write(
            "track\tsection\tmetric\t%s\n" % "\t".join(
                Stats.Summary().getHeaders()))

    # filter segments by workspace
    workspace = IO.applyIsochores(
        segments,
        annotations,
        workspaces,
        options,
        isochores,
        truncate_segments_to_workspace=options.truncate_segments_to_workspace,
        truncate_workspace_to_annotations=options.truncate_workspace_to_annotations,
        restrict_workspace=options.restrict_workspace)

    # check memory requirements
    # previous algorithm: memory requirements if all samples are stored
    # counts = segments.countsPerTrack()
    # max_counts = max(counts.values())
    # memory = 8 * 2 * options.num_samples * max_counts * len(workspace)

    # initialize sampler
    if options.sampler == "annotator":
        sampler = Engine.SamplerAnnotator(
            bucket_size=options.bucket_size,
            nbuckets=options.nbuckets)
    elif options.sampler == "shift":
        sampler = Engine.SamplerShift(
            radius=options.shift_expansion,
            extension=options.shift_extension)
    elif options.sampler == "segments":
        sampler = Engine.SamplerSegments()
    elif options.sampler == "local-permutation":
        sampler = Engine.SamplerLocalPermutation()
    elif options.sampler == "global-permutation":
        sampler = Engine.SamplerGlobalPermutation()
    elif options.sampler == "brute-force":
        sampler = Engine.SamplerBruteForce()
    elif options.sampler == "uniform":
        sampler = Engine.SamplerUniform()

    # initialize counter
    counters = []
    for counter in options.counters:
        if counter == "nucleotide-overlap":
            counters.append(Engine.CounterNucleotideOverlap())
        elif counter == "nucleotide-density":
            counters.append(Engine.CounterNucleotideDensity())
        elif counter == "segment-overlap":
            counters.append(Engine.CounterSegmentOverlap())
        elif counter == "annotation-overlap":
            counters.append(Engine.CounterAnnotationOverlap())
        elif counter == "segment-midoverlap":
            counters.append(Engine.CounterSegmentMidpointOverlap())
        elif counter == "annotation-midoverlap":
            counters.append(Engine.CounterAnnotationMidpointOverlap())
        else:
            raise ValueError("unknown counter '%s'" % counter)

    # initialize workspace generator
    if options.conditional == "unconditional":
        workspace_generator = Engine.UnconditionalWorkspace()
    elif options.conditional == "cooccurance":
        workspace_generator = Engine.ConditionalWorkspaceCooccurance()
    elif options.conditional == "annotation-centered":
        if options.conditional_expansion is None:
            raise ValueError(
                "please specify either --conditional-expansion or "
                "--conditional-extension")
        workspace_generator = Engine.ConditionalWorkspaceAnnotationCentered(
            options.conditional_extension,
            options.conditional_expansion)
    elif options.conditional == "segment-centered":
        if options.conditional_expansion is None:
            raise ValueError(
                "please specify either --conditional-expansion or "
                "--conditional-extension")

        workspace_generator = Engine.ConditionalWorkspaceSegmentCentered(
            options.conditional_extension,
            options.conditional_expansion)
    else:
        raise ValueError("unknown conditional workspace '%s'" %
                         options.conditional)

    # check if reference is compplete
    if options.reference:
        for track in segments.tracks:
            if track not in options.reference:
                raise ValueError("missing track '%s' in reference" % track)
            r = options.reference[track]
            for annotation in annotations.tracks:
                if annotation not in r:
                    raise ValueError(
                        "missing annotation '%s' in annotations for "
                        "track='%s'" % (annotation, track))

    # compute
    annotator_results = gat.run(
        segments,
        annotations,
        workspace,
        sampler,
        counters,
        workspace_generator=workspace_generator,
        num_samples=options.num_samples,
        cache=options.cache,
        outfiles=outfiles,
        output_counts_pattern=options.output_counts_pattern,
        output_samples_pattern=options.output_samples_pattern,
        sample_files=options.sample_files,
        conditional=options.conditional,
        conditional_extension=options.conditional_extension,
        reference=options.reference,
        pseudo_count=options.pseudo_count,
        num_threads=options.num_threads)

    return annotator_results


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    parser = gat.buildParser(usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    ##################################################
    description_header, descriptions, description_width = IO.readDescriptions(
        options)

    ##################################################
    size_pos, size_segment = SegmentList.getSegmentSize()
    E.debug("sizes: pos=%i segment=%i, max_coord=%i" %
            (size_pos, size_segment, 2 ** (8 * size_pos)))

    ##################################################
    # set default counter
    if not options.counters:
        options.counters.append("nucleotide-overlap")

    ##################################################
    if options.output_tables_pattern is not None:
        if "%s" not in options.output_tables_pattern:
            raise ValueError(
                "output_tables_pattern should contain at least one '%s'")

    if options.output_samples_pattern is not None:
        if "%s" not in options.output_samples_pattern:
            raise ValueError(
                "output_samples_pattern should contain at least one '%s'")

    if options.output_counts_pattern is not None:
        if "%s" not in options.output_counts_pattern:
            raise ValueError(
                "output_counts_pattern should contain at least one '%s'")

    if options.random_seed is not None:
        # initialize python random number generator
        random.seed(options.random_seed)
        # initialize numpy random number generator
        numpy.random.seed(options.random_seed)

    ##################################################
    # read fold changes that results should be compared with
    if options.null != "default":
        if not os.path.exists(options.null):
            raise OSError("file %s not found" % options.null)
        E.info("reading reference results from %s" % options.null)
        options.reference = IO.readAnnotatorResults(options.null)
    else:
        options.reference = None

    if options.input_filename_counts:
        # use pre-computed counts
        annotator_results = Engine.fromCounts(options.input_filename_counts)

    elif options.input_filename_results:
        # use previous results (re-computes fdr)
        E.info("reading gat results from %s" % options.input_filename_results)
        annotator_results = IO.readAnnotatorResults(
            options.input_filename_results)

    else:
        # do full gat analysis
        annotator_results = fromSegments(options, args)

    ##################################################
    if options.pvalue_method != "empirical":
        E.info("updating pvalues to %s" % options.pvalue_method)
        Engine.updatePValues(annotator_results, options.pvalue_method)

    ##################################################
    # output
    IO.outputResults(annotator_results,
                     options,
                     Engine.AnnotatorResultExtended.headers,
                     description_header,
                     description_width,
                     descriptions)

    IO.plotResults(annotator_results, options)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
