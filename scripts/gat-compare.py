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
gat-compare - compare two gat runs
==================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares gat runs and compares statistical significance of
fold change differences between two or more regions of interest.

This script requires counts files created by gat (option
``--output-counts-pattern``).

The reported observed value between sets 1 and 2 is the difference of
fc2 and fc1:

   * observed = fold2 - fold1

The script accepts one or more output files from a gat run.

If only a single file is given, the significance of fold change
differences are computed for each pairwise comparison of
annotations. Thus, the script might answer the question if a
transcription factor is more enriched in promotors than in enhancers.

If multiple files are given, the script computes the fold change
differences for the same annotation for all pairwise combinations of
:term:`segments of interest` across the different files. Thus, the
script might answer the question if transcription factor A is more
enriched in promotors than transcription factor B is enriched in
promotors.

Usage
-----

Example::

   python gat-compare.py tfA.counts.tsv.gz tfB.counts.tsv.gz
      
Type::

   python gat-compare.py --help

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
import math
import itertools
import numpy

import gat
import gat.Experiment as E
import gat.IOTools as IOTools
import gat.IO as IO
import gat.Engine as Engine


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-o", "--order", dest="output_order", type="choice",
                      choices=(
                          "track", "annotation", "fold", "pvalue", "qvalue", "observed"),
                      help="order results in output by fold, track, etc. [default=%default].")

    parser.add_option("-p", "--pvalue-method", dest="pvalue_method", type="choice",
                      choices=("empirical", "norm", ),
                      help="type of pvalue reported [default=%default].")

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

    parser.add_option("--pseudo-count", dest="pseudo_count", type="float",
                      help="pseudo count. The pseudo count is added to both the observed and expected overlap. "
                      " Using a pseudo-count avoids gat reporting fold changes of 0 [default=%default].")

    parser.add_option("--output-plots-pattern", dest="output_plots_pattern", type="string",
                      help="output pattern for plots [default=%default]")

    parser.set_defaults(
        pvalue_method="empirical",
        qvalue_method="BH",
        qvalue_lambda=None,
        qvalue_pi0_method="smoother",
        # pseudo count for fold change computation to avoid 0 fc
        pseudo_count=1.0,
        output_order="observed",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    input_filenames_counts = args

    ##################################################
    E.info("received %i filenames with counts" % len(input_filenames_counts))

    ##################################################
    description_header, descriptions, description_width = IO.readDescriptions(
        options)

    all_annotator_results = []

    for input_filename_counts in input_filenames_counts:

        E.info("processing %s" % input_filename_counts)

        annotator_results = gat.fromCounts(input_filename_counts)

        ##################################################
        if options.pvalue_method != "empirical":
            E.info("updating pvalues to %s" % options.pvalue_method)
            Engine.updatePValues(annotator_results, options.pvalue_method)

        ##################################################
        ##################################################
        ##################################################
        # compute global fdr
        ##################################################
        E.info("computing FDR statistics")
        Engine.updateQValues(annotator_results,
                                method=options.qvalue_method,
                                vlambda=options.qvalue_lambda,
                                pi0_method=options.qvalue_pi0_method)

        all_annotator_results.append(annotator_results)

    pseudo_count = options.pseudo_count
    results = []

    if len(all_annotator_results) == 1:
        E.info("performing pairwise comparison within a single file")

        # collect all annotations
        annotations, segments = list(), set()
        for x in all_annotator_results[0]:
            segments.add(x.track)
            annotations.append(x)

        if len(segments) != 1:
            raise NotImplementedError("multiple segments of interest")

        for data1, data2 in itertools.combinations(annotations, 2):

            # note that fold changes can be very large if there are 0 samples
            # this is fine for getting the distributional params (mean, stddev)
            fold_changes1 = data1.observed / (data1.samples + pseudo_count)
            fold_changes2 = data2.observed / (data2.samples + pseudo_count)

            # add a separate fc pseudo-count to avoid 0 values
            fold_changes1 += 0.0001
            fold_changes2 += 0.0001

            # Test is if relative fold change rfc is different from 1
            # note: rfc = fc1 / fc2 = obs1 / exp1 * obs2 / exp2
            #                       = obs1 / obs2 * exp2 / exp1
            # Thus, it is equivalent to test rfc = obs1/obs2 versus exp2 / exp1
            #
            # Convert to log space for easier plotting
            # Move the observed fold ratio in order to get an idea of the magnitude
            # of the underlying fold change
            delta_fold = data2.fold - data1.fold
            sampled_delta_fold = numpy.log(
                fold_changes1 / fold_changes2) + delta_fold
            observed_delta_fold = 0.0 + delta_fold

            result = Engine.AnnotatorResult(data1.annotation, data2.annotation,
                                               "na",
                                               observed_delta_fold,
                                               sampled_delta_fold,
                                               reference=None,
                                               pseudo_count=0)

            results.append(result)

    else:
        E.info("performing pairwise comparison between multiple files")

        ##################################################
        # perform pairwise comparison
        for index1, index2 in itertools.combinations(range(len(input_filenames_counts)), 2):
            E.info("comparing %i and %i" % (index1, index2))
            a, b = all_annotator_results[index1], all_annotator_results[index2]

            # index results in a and b
            aa = collections.defaultdict(dict)
            for x in a:
                aa[x.track][x.annotation] = x

            bb = collections.defaultdict(dict)
            for x in b:
                bb[x.track][x.annotation] = x
            
            tracks_a = set(aa.keys())
            tracks_b = set(bb.keys())
            shared_tracks = tracks_a.intersection(tracks_b)
            if len(shared_tracks) == 0:
                E.warn("no shared tracks between {} and {}".format(
                        index1, index2))
                
            for track in sorted(shared_tracks):
                E.debug("computing results for track {}".format(track))
                # get shared annotations
                annotations1 = aa[track].keys()
                annotations2 = bb[track].keys()
                shared_annotations = list(
                    set(annotations1).intersection(set(annotations2)))
                E.info("%i shared annotations" % len(shared_annotations))

                for annotation in shared_annotations:

                    # if not annotation.startswith("Ram:"): continue

                    data1 = aa[track][annotation]
                    data2 = bb[track][annotation]

                    # note that fold changes can be very large if there are 0 samples
                    # this is fine for getting the distributional params (mean,
                    # stddev)
                    fold_changes1 = data1.observed / (data1.samples + pseudo_count)
                    fold_changes2 = data2.observed / (data2.samples + pseudo_count)

                    # add a separate fc pseudo-count to avoid 0 values
                    fold_changes1 += 0.0001
                    fold_changes2 += 0.0001

                    # Test is if relative fold change rfc is different from 1
                    # note: rfc = fc1 / fc2 = obs1 / exp1 * obs2 / exp2
                    #                       = obs1 / obs2 * exp2 / exp1
                    # Thus, it is equivalent to test rfc = obs1/obs2 versus exp2 / exp1
                    #
                    # Convert to log space for easier plotting
                    # Move the observed fold ratio in order to get an idea of the magnitude
                    # of the underlying fold change
                    delta_fold = data2.fold - data1.fold
                    sampled_delta_fold = numpy.log(
                        fold_changes1 / fold_changes2) + delta_fold
                    observed_delta_fold = 0.0 + delta_fold

                    result = Engine.AnnotatorResult(track,
                                                       annotation,
                                                       "na",
                                                       observed_delta_fold,
                                                       sampled_delta_fold,
                                                       reference=None,
                                                       pseudo_count=0)

                    results.append(result)

    if len(results) == 0:
        E.critical("no results found")
        E.Stop()
        return

    IO.outputResults(results,
                     options,
                     Engine.AnnotatorResult.headers,
                     description_header,
                     description_width,
                     descriptions,
                     format_observed="%6.4f")

    IO.plotResults(results, options)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
