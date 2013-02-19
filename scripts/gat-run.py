#!/bin/env python
################################################################################
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
#################################################################################
'''
gat-run - run the genomic annotation tool
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares one or more genomic segments of interest 
against one more other genomic annotations.

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

import os, sys, re, optparse, collections, types, glob, time
import numpy

import gat
import gat.Experiment as E
import gat.IOTools as IOTools
import gat.IO as IO

def fromSegments( options, args ):
    '''run analysis from segment files. 

    This is the most common use case.
    '''

    tstart = time.time()

    segments, annotations, workspaces, isochores = IO.buildSegments( options )

    E.info( "intervals loaded in %i seconds" % (time.time() - tstart) )

    # filter segments by workspace
    workspace = IO.applyIsochores( segments, annotations, workspaces, options, isochores )

    ##################################################
    ##################################################
    ##################################################
    ## check memory requirements
    counts = segments.countsPerTrack() 
    max_counts = max(counts.values())
    # previous algorithm: memory requirements if all samples are stored
    memory = 8 * 2 * options.num_samples * max_counts * len(workspace)

    ##################################################
    ##################################################
    ##################################################
    # initialize sampler
    if options.sampler == "annotator":
        sampler = gat.SamplerAnnotator(
            bucket_size = options.bucket_size,
            nbuckets = options.nbuckets )
    elif options.sampler == "shift":
        sampler = gat.SamplerShift( 
            radius = options.shift_expansion,
            extension = options.shift_extension )
    elif options.sampler == "segments":
        sampler = gat.SamplerSegments()
    elif options.sampler == "local-permutation":
        sampler = gat.SamplerLocalPermutation()
    elif options.sampler == "global-permutation":
        sampler = gat.SamplerGlobalPermutation()
    elif options.sampler == "brute-force":
        sampler = gat.SamplerBruteForce()
    elif options.sampler == "uniform":
        sampler = gat.SamplerUniform()
        
    ##################################################
    ##################################################
    ##################################################
    # initialize counter
    counters = []
    for counter in options.counters:
        if counter == "nucleotide-overlap":
            counters.append( gat.CounterNucleotideOverlap() )
        elif counter == "nucleotide-density":
            counters.append( gat.CounterNucleotideDensity() )
        elif counter == "segment-overlap":
            counters.append( gat.CounterSegmentOverlap() )
        elif counter == "annotations-overlap":
            counters.append( gat.CounterAnnotationsOverlap() )
        elif counter == "segment-midoverlap":
            counters.append( gat.CounterSegmentMidpointOverlap() )
        elif counter == "annotations-midoverlap":
            counters.append( gat.CounterAnnotationsMidpointOverlap() )
        else:
            raise ValueError("unknown counter '%s'" % counter )

    ##################################################
    ##################################################
    ##################################################
    ## initialize workspace generator
    if options.conditional == "unconditional":
        workspace_generator = gat.UnconditionalWorkspace()
    elif options.conditional == "cooccurance":
        workspace_generator = gat.ConditionalWorkspaceCooccurance()
    elif options.conditional == "annotation-centered":
        workspace_generator = gat.ConditionalWorkspaceAnnotationCentered( options.conditional_extension,
                                                                          options.conditional_expansion )
    elif options.conditional == "segment-centered":
        workspace_generator = gat.ConditionalWorkspaceSegmentCentered( options.conditional_extension,
                                                                       options.conditional_expansion )
    else:
        raise ValueError("unknown conditional workspace '%s'" % options.conditional )


    ##################################################
    ##################################################
    ##################################################
    ## check if reference is compplete
    ##################################################
    if options.reference:
        reference = options.reference
        for track in segments.tracks:
            if track not in options.reference:
                raise ValueError("missing track '%s' in reference" % track )
            r = options.reference[track]
            for annotation in annotations.tracks:
                if annotation not in r:
                    raise ValueError("missing annotation '%s' in annotations for track='%s'" % (annotation, track ))

    ##################################################
    ##################################################
    ##################################################
    ## open sample stats outfile
    ##################################################
    if "sample" in options.output_stats or \
            "all" in options.output_stats or \
                len( [ x for x in options.output_stats if re.search( x, "sample" ) ] ) > 0:
        outfile_sample_stats = E.openOutputFile( "sample_stats" )
    else:
        outfile_sample_stats = None

    

    ##################################################
    ##################################################
    ##################################################
    ## compute
    ##################################################
    annotator_results = gat.run( segments, 
                                 annotations, 
                                 workspace,
                                 sampler, 
                                 counters,
                                 workspace_generator = workspace_generator,
                                 num_samples = options.num_samples,
                                 cache = options.cache,
                                 outfile_sample_stats = outfile_sample_stats,
                                 output_counts_pattern = options.output_counts_pattern,
                                 output_samples_pattern = options.output_samples_pattern,
                                 sample_files = options.sample_files,
                                 conditional = options.conditional,
                                 conditional_extension = options.conditional_extension,
                                 reference = options.reference,
                                 pseudo_count = options.pseudo_count,
                                 num_threads = options.num_threads )

    return annotator_results

def fromResults( filename ):
    '''load annotator results from a tab-separated results table.'''

    annotator_results = collections.defaultdict( dict )

    with IOTools.openFile(filename, "r") as infile:
        for line in infile:
            if line.startswith("#"): continue
            if line.startswith("track"): continue
            r = gat.DummyAnnotatorResult._fromLine( line ) 
            annotator_results[r.track][r.annotation] = r
            
    return annotator_results

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--annotation-file", "--annotations", dest="annotation_files", type="string", action="append",
                      help="filename with annotations [default=%default]."  )

    parser.add_option("-s", "--segment-file", "--segments", dest="segment_files", type="string", action="append",
                      help="filename with segments. Also accepts a glob in parentheses [default=%default]."  )

    parser.add_option("-w", "--workspace-file", "--workspace", dest="workspace_files", type="string", action="append",
                      help="filename with workspace segments. Also accepts a glob in parentheses [default=%default]."  )

    parser.add_option("-i", "--isochore-file", "--isochores", dest="isochore_files", type="string", action="append",
                      help="filename with isochore segments. Also accepts a glob in parentheses [default=%default]."  )

    parser.add_option("-l", "--sample-file", dest="sample_files", type="string", action="append",
                      help="filename with sample files. Start processing from samples [default=%default]."  )

    parser.add_option("-c", "--counter", dest="counters", type="choice", action="append",
                      choices=("nucleotide-overlap", 
                               "nucleotide-density",
                               "segment-overlap", 
                               "segment-midoverlap",
                               "annotations-overlap", 
                               "annotations-midoverlap" ),
                      help="quantity to use for estimating enrichment [default=%default]."  )

    parser.add_option("-m", "--sampler", dest="sampler", type="choice",
                      choices=("annotator", 
                               "segments",
                               "shift",
                               "local-permutation",
                               "global-permutation",
                               "uniform",
                               "brute-force"),
                      help="quantity to test [default=%default]."  )

    parser.add_option("-n", "--num-samples", dest="num_samples", type="int", 
                      help="number of samples to compute [default=%default]."  )

    parser.add_option("-e", "--cache", dest="cache", type="string", 
                      help="filename for caching samples [default=%default]."  )

    parser.add_option("-o", "--order", dest="output_order", type="choice",
                      choices = ( "track", "annotation", "fold", "pvalue", "qvalue" ),
                      help="order results in output by fold, track, etc. [default=%default]."  )

    parser.add_option("-p", "--pvalue-method", dest="pvalue_method", type="choice",
                      choices = ( "empirical", "norm", ),
                      help="type of pvalue reported [default=%default]."  )

    parser.add_option("-q", "--qvalue-method", dest="qvalue_method", type="choice",
                      choices = ( "storey", "BH", "bonferroni", "holm", "hommel", "hochberg", "BY", "none" ),
                      help="method to perform multiple testing correction by controlling the fdr [default=%default]."  )

    parser.add_option("-t", "--num-threads", dest="num_threads", type="int",
                      help = "number of threads to use for sampling [default=%default]" )

    parser.add_option( "--qvalue-lambda", dest="qvalue_lambda", type="float",
                      help="fdr computation: lambda [default=%default]."  )

    parser.add_option( "--qvalue-pi0-method", dest="qvalue_pi0_method", type="choice",
                       choices = ("smoother", "bootstrap" ),
                       help="fdr computation: method for estimating pi0 [default=%default]."  )

    parser.add_option( "--input-counts-file", dest="input_filename_counts", type="string", 
                      help="start processing from counts - no segments required [default=%default]."  )

    parser.add_option( "--input-results-file", dest="input_filename_results", type="string", 
                      help="start processing from results - no segments required [default=%default]."  )

    parser.add_option( "--output-tables-pattern", dest="output_tables_pattern", type="string", 
                      help="output pattern for result tables. Used if there are multiple counters used [default=%default]."  )

    parser.add_option( "--output-counts-pattern", dest="output_counts_pattern", type="string", 
                      help="output pattern for counts [default=%default]."  )

    parser.add_option( "--output-plots-pattern", dest="output_plots_pattern", type="string", 
                       help="output pattern for plots [default=%default]" )

    parser.add_option( "--output-samples-pattern", dest="output_samples_pattern", type="string", 
                       help="output pattern for samples. Samples are stored in bed format, one for "
                            " each segment [default=%default]" )

    parser.add_option( "--descriptions", dest="input_filename_descriptions", type="string", 
                       help="filename mapping annotation terms to descriptions. "
                            " if given, the output table will contain additional columns "
                            " [default=%default]" )

    parser.add_option( "--bucket-size", dest="bucket_size", type="int", 
                       help="size of a bin for histogram of segment lengths [default=%default]" )

    parser.add_option( "--nbuckets", dest="nbuckets", type="int", 
                       help="number of bins for histogram of segment lengths [default=%default]" )

    parser.add_option( "--output-stats", dest="output_stats", type="choice", action="append",
                       choices = ( "all", 
                                   "annotations", "segments", 
                                   "workspaces", "isochores",
                                   "overlap" ),
                       help="output overlap summary stats [default=%default]."  )

    parser.add_option( "--output-bed", dest="output_bed", type="choice", action="append",
                       choices = ( "all", 
                                   "annotations", "segments", 
                                   "workspaces", "isochores",
                                   "overlap" ),
                       help="output bed files [default=%default]."  )
    
    parser.add_option( "--ignore-segment-tracks", dest="ignore_segment_tracks", action="store_true", 
                       help="ignore segment tracks - all segments belong to one track [default=%default]" )

    parser.add_option( "--enable-split-tracks", dest="enable_split_tracks", action="store_true", 
                       help="permit the same track to be in multiple files [default=%default]" )

    parser.add_option( "--conditional", dest="conditional", type="choice",
                       choices = ( "unconditional", "annotation-centered", "segment-centered", "cooccurance" ),
                       help="conditional workspace creation [default=%default]"
                       " cooccurance - compute enrichment only within workspace segments that contain both segments "
                       " and annotations"
                       " annotation-centered - workspace centered around annotations. See --conditional-extension"
                       " segment-centered - workspace centered around segments. See --conditional-extension" )

    parser.add_option( "--conditional-extension", dest="conditional_extension", type="int",
                      help="if workspace is created conditional, extend by this amount (in bp) [default=%default]."  )

    parser.add_option( "--conditional-expansion", dest="conditional_expansion", type="float",
                      help="if workspace is created conditional, expand by this amount (ratio) [default=%default]."  )

    parser.add_option( "--restrict-workspace", dest="restrict_workspace", action="store_true", 
                       help="restrict workspace to those segments that contain both track"
                       " and annotations [default=%default]" )

    parser.add_option( "--truncate-workspace-to-annotations", dest="truncate_workspace_to_annotations", action="store_true", 
                       help="truncate workspace with annotations [default=%default]" )

    parser.add_option( "--shift-extension", dest="shift_extension", type="float",
                      help="if the sampling method is 'shift', create a segment of size # anound the segment"
                       " to determine the size of the region for shifthing [default=%default]."  )

    parser.add_option( "--shift-expansion", dest="shift_expansion", type="float",
                      help="if the sampling method is 'shift', multiply each segment by # "
                           " to determine the size of the region for shifthing [default=%default]."  )

    parser.add_option( "--pseudo-count", dest="pseudo_count", type="float",
                      help="pseudo count. The pseudo count is added to both the observed and expected overlap. "
                       " Using a pseudo-count avoids gat reporting fold changes of 0 [default=%default]."  )

    parser.add_option( "--null", dest="null", type="string",
                      help="null hypothesis. The default is to test categories for enrichment/depletion. "
                           " If a filename with gat output is given, gat will test for the difference "
                           " in fold change between the segments supplied and in the other file [default=%default]."  )

    parser.set_defaults(
        annotation_files = [],
        segment_files = [],
        workspace_files = [],  
        sample_files = [],
        num_samples = 1000,
        nbuckets = 100000,
        bucket_size = 1,
        counters = [],
        output_stats = [],
        output_bed = [],
        output_order = "fold",
        cache = None,
        input_filename_counts = None,
        input_filename_results = None,
        pvalue_method = "empirical",
        output_tables_pattern = "%s.tsv.gz",
        output_counts_pattern = None,
        output_plots_pattern = None,
        output_samples_pattern = None,
        qvalue_method = "BH",
        qvalue_lambda = None,
        qvalue_pi0_method = "smoother",
        sampler = "annotator",
        ignore_segment_tracks = False,
        input_filename_descriptions = None,
        conditional = "unconditional",
        conditional_extension = None,
        conditional_expansion = None,
        restrict_workspace = False,
        truncate_workspace_to_annotations = False,
        enable_split_tracks = False,
        shift_expansion = 2.0,
        shift_extension = 0,
        # pseudo count for fold change computation to avoid 0 fc
        pseudo_count = 1.0,
        null = "default",
        num_threads = 1,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ##################################################
    description_header, descriptions, description_width = IO.readDescriptions( options )

    ##################################################
    size_pos, size_segment = gat.csegmentlist.getSegmentSize()
    E.debug( "sizes: pos=%i segment=%i, max_coord=%i" % (size_pos, size_segment, 2**(8 * size_pos )))

    ##################################################
    # set default counter
    if not options.counters:
        options.counters.append( "nucleotide-overlap" )

    ##################################################
    if options.output_tables_pattern != None: 
        if "%s" not in options.output_tables_pattern:
            raise ValueError( "output_tables_pattern should contain at least one '%s'")

    if options.output_samples_pattern != None: 
        if "%s" not in options.output_samples_pattern:
            raise ValueError( "output_samples_pattern should contain at least one '%s'")

    if options.output_counts_pattern != None: 
        if "%s" not in options.output_counts_pattern:
            raise ValueError( "output_counts_pattern should contain at least one '%s'")


    ##################################################
    # read fold changes that results should be compared with
    if options.null != "default":
        if not os.path.exists( options.null ):
            raise OSError( "file %s not found" % options.null )
        E.info( "reading reference results from %s" % options.null )
        options.reference = fromResults( options.null )
    else:
        options.reference = None

    if options.input_filename_counts:
        # use pre-computed counts
        annotator_results = gat.fromCounts( options.input_filename_counts )

    elif options.input_filename_results:
        # use previous results (re-computes fdr)
        E.info( "reading gat results from %s" % options.input_filename_results )
        annotator_results = fromResults( options.input_filename_results )

    else:
        # do full gat analysis
        annotator_results = fromSegments( options, args )

    ##################################################
    if options.pvalue_method != "empirical":
        E.info("updating pvalues to %s" % options.pvalue_method )
        gat.updatePValues( annotator_results, options.pvalue_method )

    ##################################################
    ## output
    IO.outputResults( annotator_results, 
                      options, 
                      gat.AnnotatorResultExtended.headers,
                      description_header, 
                      description_width,
                      descriptions )

    IO.plotResults( annotator_results, options )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
