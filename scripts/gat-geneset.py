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




Usage
-----

Example::

   python gatrun.py 
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
import csegmentlist

import scipy.stats
import Bed

from bx.intervals.intersection import Intersecter, Interval

GENESET_RESULT = collections.namedtuple("GENESET_RESULT",
                                        ("track",
                                         "annotation",
                                         "genes_in_workspace",
                                         "genes_with_annotation",
                                         "genes_with_segments",
                                         "genes_with_annotation_and_segments",
                                         "expected",
                                         "fold",
                                         "pvalue",
                                         "qvalue"))


def main(argv):

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-a", "--gene-file", "--annotations", dest="annotation_files", type="string", action="append",
                      help="filename with annotations - here, location of genes [default=%default].")

    parser.add_option("-s", "--segment-file", "--segments", dest="segment_files", type="string", action="append",
                      help="filename with segments. Also accepts a glob in parentheses [default=%default].")

    parser.add_option("-w", "--workspace-file", "--workspace", dest="workspace_files", type="string", action="append",
                      help="filename with workspace segments. Also accepts a glob in parentheses [default=%default].")

    parser.add_option("-g", "--number-of-genes", dest="number_of_genes", type="int",
                      help="total number of genes [default=%default]")

    parser.add_option("-m", "--annotation-file", dest="annotation_file", type="string",
                      help="filename mapping genes to annotations [default=%default]")

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

    parser.set_defaults(
        annotation_files=[],
        segment_files=[],
        workspace_files=[],
        sample_files=[],
        annotation_file=None,
        num_samples=1000,
        nbuckets=100000,
        bucket_size=1,
        counter="nucleotide-overlap",
        output_stats=[],
        output_bed=[],
        output_filename_counts=None,
        output_order="fold",
        cache=None,
        input_filename_counts=None,
        input_filename_results=None,
        pvalue_method="empirical",
        output_plots_pattern=None,
        output_samples_pattern=None,
        qvalue_method="storey",
        qvalue_lambda=None,
        qvalue_pi0_method="smoother",
        sampler="annotator",
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
        number_of_genes=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    tstart = time.time()

    # load segments
    options.segment_files = IO.expandGlobs(options.segment_files)
    options.annotation_files = IO.expandGlobs(options.annotation_files)
    options.workspace_files = IO.expandGlobs(options.workspace_files)

    # read one or more segment files
    segments = IO.readSegmentList("segments", options.segment_files, options)
    if options.ignore_segment_tracks:
        segments.merge(delete=True)
        E.info("merged all segments into one track with %i segments" %
               len(segments))

    if len(segments) > 1000:
        raise ValueError(
            "too many (%i) segment files - use track definitions or --ignore-segment-tracks" % len(segments))

    # load workspace
    workspaces = IO.readSegmentList(
        "workspaces", options.workspace_files, options, options.enable_split_tracks)

    # intersect workspaces to build a single workspace
    E.info("collapsing workspaces")
    workspaces.collapse()

    # use merged workspace only, discard others
    workspaces.restrict("collapsed")
    workspace = workspaces["collapsed"]

    E.info("intervals loaded in %i seconds" % (time.time() - tstart))

    ############################################
    # load table mapping a gene id to annotations
    gene2annotations = IOTools.readMultiMap(IOTools.openFile(options.annotation_file),
                                            has_header=True)
    annotations = set([y for x in gene2annotations.values() for y in x])
    E.info("loaded %i annotations for %i genes" %
           (len(gene2annotations), len(annotations)))

    ############################################
    # load bed file with gene coordinates
    assert len(options.annotation_files) == 1
    indexed_genes = collections.defaultdict(Intersecter)
    total_genes = 0
    # number of genes per contig
    contig2ngenes = collections.defaultdict(int)
    # compute number of genes with a particular annotation
    # per contig
    annotation2ngenes = collections.defaultdict(int)
    for line in IOTools.openFile(options.annotation_files[0]):
        if line.startswith("#"):
            continue
        contig, start, end, gene_id = line[:-1].split("\t")[:4]
        indexed_genes[contig].add_interval(
            Interval(int(start), int(end), gene_id))
        contig2ngenes[contig] += 1
        total_genes += 1
        try:
            for annotation in gene2annotations[gene_id]:
                annotation2ngenes[annotation] += 1
        except KeyError:
            pass
    E.info("indexed locations for %i contigs" % len(indexed_genes))

    ############################################
    description_header, descriptions, description_width = IO.readDescriptions(
        options)

    ############################################
    ############################################
    # compute results
    E.info("computing counts")

    results = []
    # iterate over segments
    for segment, segmentdict in segments.iteritems():

        # genes hit by segments per annotation
        genes_hit_by_segments_with_annotations = collections.defaultdict(int)

        # genes hit by segments
        genes_hit_by_segments = 0

        for contig, ss in segmentdict.iteritems():
            for start, end in ss:
                overlapping_genes = list(
                    indexed_genes[contig].find(start, end))
                genes_hit_by_segments += len(overlapping_genes)
                for x in overlapping_genes:
                    gene_id = x.value
                    try:
                        for annotation in gene2annotations[gene_id]:
                            genes_hit_by_segments_with_annotations[
                                annotation] += 1
                    except KeyError:
                        pass

        # N = number of genes in genome
        N = total_genes
        # n   = number of genes selected by segments
        n = genes_hit_by_segments

        for annotation in annotations:
            # K = number of genes carrying annotation
            K = annotation2ngenes[annotation]
            # k = number of genes selected by segments and with annotation
            k = genes_hit_by_segments_with_annotations[annotation]

            if n == 0 or N == 0 or K == 0:
                expected = 0
                fold = 1.0
                pvalue = 1.0
            else:
                expected = float(n * K) / N
                fold = k / expected
                pvalue = scipy.stats.hypergeom.sf(k - 1, N, K, n)

            r = GENESET_RESULT._make((
                segment, annotation,
                N,
                K,
                n,
                k,
                expected,
                fold,
                pvalue,
                1.0))

            results.append(r)

    IO.outputResults(results,
                     options,
                     GENESET_RESULT._fields,
                     description_header,
                     description_width,
                     descriptions)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
