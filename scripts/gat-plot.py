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
gat-plot - plot results from a gat analysis
===========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes the results of a ``gat-run.py` or ``gat-compare.py``
and plots the results.

This script requires matplotlib.

Usage
-----

Example::

   python gat-plot.py --input-filename-results=gat.results.tsv.gz
   python gat-plot.py --input-filename-counts=gat.counts.tsv.gz
 
Type::

   python gatplot.py --help

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

import gat
import gat.Experiment as E
import gat.IOTools as IOTools
import gat.IO as IO

try:
    import matplotlib.pyplot as plt
    HASPLOT = True
except (ImportError, RuntimeError):
    HASPLOT = False


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
        x.observed, x.expected, x.lower95, x.upper95, x.stddev, x.fold, x.pvalue, x.qvalue = \
            map(float, data[2:])
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


def buildPlotFilename(options, key):
    filename = re.sub("%s", key, options.output_plots_pattern)
    filename = re.sub("[^a-zA-Z0-9-_./]", "_", filename)
    dirname = os.path.dirname(filename)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
    return filename


def plotBarplots(annotator_results, options):
    '''output a series of bar-plots.

    Output for each track.

    Significant results are opaque, while
    non-significant results are transparent.'''

    for track in annotator_results:
        plt.figure()
        r = annotator_results[track]
        keys, values = zip(*r.items())
        pos = range(len(r))
        bars = plt.barh(pos, [x.fold for x in values])
        for b, v in zip(bars, values):
            if v.qvalue > 0.05:
                b.set_alpha(0.10)

        filename = buildPlotFilename(options, "bars-%s" % track)
        plt.yticks(pos, keys)
        plt.axvline(x=1, color="r")
        plt.savefig(filename)


def plotBarplot(annotator_results, options):
    '''output a single bar-plots.

    Output for each track.

    Significant results are opaque, while
    non-significant results are transparent.'''

    ntracks = len(annotator_results)
    height = 1.0 / float(ntracks)

    plt.figure()

    for trackid, track in enumerate(annotator_results):

        r = annotator_results[track]
        rr = r.items()
        rr.sort()
        keys, values = zip(*rr)
        pos = numpy.arange(0, len(r), 1) + trackid * height
        bars = plt.barh(pos,
                        [x.fold for x in values],
                        height=height,
                        label=track,
                        xerr=[x.stddev / x.expected for x in values],
                        color="bryg"[trackid % 4])
        for b, v in zip(bars, values):
            if v.pvalue > 0.05:
                b.set_alpha(0.10)

    pos = range(len(r))

    plt.yticks(pos, keys)
    plt.axvline(x=1, color="r")
    filename = buildPlotFilename(options, "bars-all")
    plt.legend()
    plt.savefig(filename)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-l", "--sample-file", dest="sample_files", type="string", action="append",
                      help="filename with sample files. Start processing from samples [default=%default].")

    parser.add_option("-o", "--order", dest="output_order", type="choice",
                      choices=(
                          "track", "annotation", "fold", "pvalue", "qvalue"),
                      help="order results in output by fold, track, etc. [default=%default].")

    parser.add_option("-p", "--pvalue-method", dest="pvalue_method", type="choice",
                      choices=("empirical", "norm", ),
                      help="type of pvalue reported [default=%default].")

    parser.add_option("--results-file", dest="input_filename_results", type="string",
                      help="start processing from results - no segments required [default=%default].")

    parser.add_option("--output-plots-pattern", dest="output_plots_pattern", type="string",
                      help="output pattern for plots [default=%default]")

    parser.add_option("--output-samples-pattern", dest="output_samples_pattern", type="string",
                      help="output pattern for samples. Samples are stored in bed format, one for "
                      " each segment [default=%default]")

    parser.add_option("--plots", dest="plots", type="choice",
                      choices=("all",
                               "bars-per-track",
                               "bars", ),
                      help="plots to be created [default=%default].")

    parser.set_defaults(
        sample_files=[],
        num_samples=1000,
        output_stats=[],
        output_filename_counts=None,
        output_order="fold",
        input_filename_results=None,
        pvalue_method="empirical",
        output_plots_pattern=None,
        plots=[],
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    annotator_results = IO.readAnnotatorResults(options.input_filename_results)

    if "speparate-bars" in options.plots:
        plotBarplots(annotator_results, options)
    if "bars" in options.plots:
        plotBarplot(annotator_results, options)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
