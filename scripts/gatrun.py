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

import os, sys, re, optparse, collections, types, glob, time
import numpy

import gat
import gat.Experiment as E
import gat.IOTools as IOTools

try:
    import matplotlib.pyplot as plt
    HASPLOT = True
except (ImportError,RuntimeError):
    HASPLOT = False

def fromSegments( options, args ):
    '''run analysis from segment files. 

    This is the most common use case.
    '''

    tstart = time.time()

    def expandGlobs( infiles ):
        return IOTools.flatten( [ glob.glob( x ) for x in infiles ] )
        
    options.segment_files = expandGlobs( options.segment_files )
    options.annotation_files = expandGlobs( options.annotation_files )
    options.workspace_files = expandGlobs( options.workspace_files )
    options.sample_files = expandGlobs( options.sample_files )

    ##################################################
    # arguments sanity check
    if not options.segment_files:
        raise ValueError("please specify at least one segment file" )
    if not options.annotation_files:
        raise ValueError("please specify at least one annotation file" )
    if not options.workspace_files:
        raise ValueError("please specify at least one workspace file" )

    ##################################################
    ##################################################
    ##################################################
    # process input
    def dumpStats( coll, section ):
        if section in options.output_stats or \
                "all" in options.output_stats or \
                len( [ x for x in options.output_stats if re.search( x, section ) ] ) > 0:
            coll.outputStats( E.openOutputFile( section ) )

    def dumpBed( coll, section ):
        if section in options.output_bed or \
                "all" in options.output_bed or \
                len( [ x for x in options.output_bed if re.search( x, section ) ] ) > 0:
            coll.save( E.openOutputFile( section + ".bed" ) )

    def readSegmentList( label, filenames, enable_split_tracks = False ):
        # read one or more segment files
        results = gat.IntervalCollection( name = label )
        E.info( "%s: reading tracks from %i files" % (label, len(filenames)))
        results.load( filenames, split_tracks = enable_split_tracks )
        E.info( "%s: read %i tracks from %i files" % (label, len(results), len(filenames)))
        dumpStats( results, "stats_%s_raw" % label )
        results.normalize()
        dumpStats( results, "stats_%s_normed" % label )
        return results

    # read one or more segment files
    segments = readSegmentList( "segments", options.segment_files)
    if options.ignore_segment_tracks:
        segments.merge( delete = True)
        E.info( "merged all segments into one track with %i segments" % len(segments))

    if len(segments) > 1000: 
        raise ValueError( "too many (%i) segment files - use track definitions or --ignore-segment-tracks" % len(segments) )
    
    annotations = readSegmentList( "annotations", options.annotation_files, options.enable_split_tracks )
    workspaces = readSegmentList( "workspaces", options.workspace_files, options.enable_split_tracks )

    # intersect workspaces to build a single workspace
    E.info( "collapsing workspaces" )
    dumpStats( workspaces, "stats_workspaces_input" )
    workspaces.collapse()
    dumpStats( workspaces, "stats_workspaces_collapsed" )

    # use merged workspace only, discard others
    workspaces.restrict("collapsed")

    # build isochores or intersect annotations/segments with workspace
    if options.isochore_files:
        
        # read one or more isochore files
        isochores = gat.IntervalCollection( name = "isochores" )
        E.info( "%s: reading isochores from %i files" % ("isochores", len(options.isochore_files)))
        isochores.load( options.isochore_files )
        dumpStats( isochores, "stats_isochores_raw" )

        # merge isochores and check if consistent (fully normalized)
        isochores.sort()

        # check that there are no overlapping segments within isochores
        isochores.check()

        # TODO: flag is_normalized not properly set
        isochores.normalize()

        # check that there are no overlapping segments between isochores

        # truncate isochores to workspace
        # crucial if isochores are larger than workspace.
        isochores.intersect( workspaces["collapsed"] )

        # intersect isochores and workspaces, segments and annotations
        E.info( "adding isochores to workspace" )
        workspaces.toIsochores( isochores )
        annotations.toIsochores( isochores )
        segments.toIsochores( isochores )
        
        if workspaces.sum() == 0:
            raise ValueError( "isochores and workspaces do not overlap" )
        if annotations.sum() == 0:
            raise ValueError( "isochores and annotations do not overlap" )
        if segments.sum() == 0:
            raise ValueError( "isochores and segments do not overlap" )

        dumpStats( workspaces, "stats_workspaces_isochores" )
        dumpStats( annotations, "stats_annotations_isochores" )
        dumpStats( segments, "stats_segments_isochores" )
    
        dumpBed( workspaces, "workspaces_isochores" )
        dumpBed( annotations, "annotations_isochores" )
        dumpBed( segments, "segments_isochores" )

    else:
        # intersect workspace and segments/annotations
        annotations.filter( workspaces["collapsed"] )
        segments.filter( workspaces["collapsed"] )
        
        dumpStats( annotations, "stats_annotations_pruned" )
        dumpStats( segments, "stats_segments_pruned" )

    workspace = workspaces["collapsed"] 

    if options.restrict_workspace:
        E.info( "restricting workspace" )
        # this is very cumbersome - refactor merge and collapse
        # to return an IntervalDictionary instead of adding it
        # to the list of tracks
        for x in (segments, annotations):
            if "merged" in segments:
                workspace.filter( segments["merged"] )
            else:
                segments.merge()
                workspace.filter( segments["merged"] )
                del segments[merged]

        dumpStats( workspaces, "stats_workspaces_restricted" )
        
    # segments.dump( open("segments_dump.bed", "w" ) )
    # workspaces.dump( open("workspaces_dump.bed", "w" ) )

    E.info( "intervals loaded in %i seconds" % (time.time() - tstart) )

    # output overlap stats
    # output segment densities per workspace
    for track in segments.tracks:
        workspaces.outputOverlapStats( E.openOutputFile( "overlap_%s" % track), 
                                       segments[track] )

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
        
    ##################################################
    ##################################################
    ##################################################
    # initialize counter
    if options.counter == "nucleotide-overlap":
        counter = gat.CounterNucleotideOverlap()
    elif options.counter == "nucleotide-density":
        counter = gat.CounterNucleotideDensity()
    elif options.counter == "segment-overlap":
        counter = gat.CounterSegmentOverlap()
    else:
        raise ValueError("unknown counter '%s'" % options.counter )

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
                                 counter,
                                 workspace_generator = workspace_generator,
                                 num_samples = options.num_samples,
                                 cache = options.cache,
                                 outfile_sample_stats = outfile_sample_stats,
                                 output_counts = options.output_filename_counts,
                                 output_samples_pattern = options.output_samples_pattern,
                                 sample_files = options.sample_files,
                                 conditional = options.conditional,
                                 conditional_extension = options.conditional_extension )

    return annotator_results

def fromResults( filename ):
    '''load annotator results from a tab-separated results table.'''

    annotator_results = collections.defaultdict( dict )

    with open(filename, "r") as infile:
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

    parser.add_option("-c", "--counter", dest="counter", type="choice",
                      choices=("nucleotide-overlap", 
                               "nucleotide-density",
                               "segment-overlap", ),
                      help="quantity to test [default=%default]."  )

    parser.add_option("-m", "--sampler", dest="sampler", type="choice",
                      choices=("annotator", 
                               "segments",
                               "shift" ),
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
                      choices = ( "storey", ),
                      help="method to perform multiple testing correction by controlling the fdr [default=%default]."  )

    parser.add_option( "--qvalue-lambda", dest="qvalue_lambda", type="float",
                      help="fdr computation: lambda [default=%default]."  )

    parser.add_option( "--qvalue-pi0-method", dest="qvalue_pi0_method", type="choice",
                       choices = ("smoother", "bootstrap" ),
                       help="fdr computation: method for estimating pi0 [default=%default]."  )

    parser.add_option( "--counts-file", dest="input_filename_counts", type="string", 
                      help="start processing from counts - no segments required [default=%default]."  )

    parser.add_option( "--output-counts-file", dest="output_filename_counts", type="string", 
                      help="output counts to filename [default=%default]."  )

    parser.add_option( "--results-file", dest="input_filename_results", type="string", 
                      help="start processing from results - no segments required [default=%default]."  )

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

    parser.add_option( "--shift-extension", dest="shift_extension", type="float",
                      help="if the sampling method is 'shift', create a segment of size # anound the segment"
                       " to determine the size of the region for shifthing [default=%default]."  )

    parser.add_option( "--shift-expansion", dest="shift_expansion", type="float",
                      help="if the sampling method is 'shift', multiply each segment by # "
                           " to determine the size of the region for shifthing [default=%default]."  )

    parser.set_defaults(
        annotation_files = [],
        segment_files = [],
        workspace_files = [],  
        sample_files = [],
        num_samples = 1000,
        nbuckets = 100000,
        bucket_size = 1,
        counter = "nucleotide-overlap",
        output_stats = [],
        output_bed = [],
        output_filename_counts = None,
        output_order = "fold",
        cache = None,
        input_filename_counts = None,
        input_filename_results = None,
        pvalue_method = "empirical",
        output_plots_pattern = None,
        output_samples_pattern = None,
        qvalue_method = "storey",
        qvalue_lambda = None,
        qvalue_pi0_method = "smoother",
        sampler = "annotator",
        ignore_segment_tracks = False,
        input_filename_descriptions = None,
        conditional = "unconditional",
        conditional_extension = None,
        conditional_expansion = None,
        restrict_workspace = False,
        enable_split_tracks = False,
        shift_expansion = 2.0,
        shift_extension = 0,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ##################################################
    description_header, descriptions, description_width = [], {}, 0
    if options.input_filename_descriptions:
        E.info( "reading descriptions from %s" % options.input_filename_descriptions )

        with IOTools.openFile( options.input_filename_descriptions ) as inf:
            first = True
            for line in inf:
                if line.startswith("#"): continue
                data = line[:-1].split( "\t" )

                if description_width: assert len(data) -1 == description_width
                else: description_width = len(data) - 1

                if first: 
                    description_header = data[1:]
                    first = False
                else:
                    descriptions[data[0]] = data[1:]

    size_pos, size_segment = gat.getSegmentSize()
    E.debug( "sizes: pos=%i segment=%i, max_coord=%i" % (size_pos, size_segment, 2**(8 * size_pos )))

    ##################################################
    if options.input_filename_counts:
        annotator_results = gat.fromCounts( options.input_filename_counts )
    elif options.input_filename_results:
        E.info( "reading annotator results from %s" % options.input_filename_results )
        annotator_results = fromResults( options.input_filename_results )
    else:
        annotator_results = fromSegments( options, args )

    if options.pvalue_method != "empirical":
        E.info("updating pvalues to %s" % options.pvalue_method )
        gat.updatePValues( gat.iterator_results(annotator_results), options.pvalue_method )

    ##################################################
    ##################################################
    ##################################################
    ## compute global fdr
    ##################################################
    E.info( "computing FDR statistics" )
    gat.updateQValues( list(gat.iterator_results(annotator_results)), 
                       method = options.qvalue_method,
                       vlambda = options.qvalue_lambda,
                       pi0_method = options.qvalue_pi0_method )

    ##################################################
    ##################################################
    ## output
    ##################################################
    outfile = options.stdout

    outfile.write( "\t".join( gat.AnnotatorResult.headers + description_header ) + "\n" )

    output = list( gat.iterator_results( annotator_results ) )
    if options.output_order == "track":
        output.sort( key = lambda x: (x.track, x.annotation) )
    elif options.output_order == "annotation":
        output.sort( key = lambda x: (x.annotation, x.track) )
    elif options.output_order == "fold":
        output.sort( key = lambda x: x.fold )
    elif options.output_order == "pvalue":
        output.sort( key = lambda x: x.pvalue )
    elif options.output_order == "qvalue":
        output.sort( key = lambda x: x.qvalue )
    else:
        raise ValueError("unknown sort order %s" % options.output_order )

    for result in output:
        outfile.write( str(result) )
        if descriptions:
            try:
                outfile.write( "\t" + "\t".join( descriptions[result.annotation] ) )
            except KeyError:
                outfile.write( "\t" + "\t".join( [""] * description_width ) )
        outfile.write("\n")
    
    ##################################################
    # plot histograms
    if options.output_plots_pattern and HASPLOT:

        def buildPlotFilename( options, key ):
            filename = re.sub("%s", key, options.output_plots_pattern)
            filename = re.sub("[^a-zA-Z0-9-_./]", "_", filename )
            dirname = os.path.dirname( filename )
            if dirname and not os.path.exists( dirname ): os.makedirs( dirname )
            return filename

        E.info("plotting sample stats" )

        for r in gat.iterator_results(annotator_results):
            plt.figure()
            key = "%s-%s" % (r.track, r.annotation)
            s = r.samples
            hist, bins = numpy.histogram( s,
                                          bins = 100 )
            
            # convert to density
            hist = numpy.array( hist, dtype = numpy.float )
            hist /= sum(hist)

            # plot bars
            plt.bar( bins[:-1], hist, width=1.0, label = key )
            
            # plot estimated 
            sigma = r.stddev
            mu = r.expected
            plt.plot(bins, 
                     1.0/(sigma * numpy.sqrt(2 * numpy.pi)) *
                     numpy.exp( - (bins - mu)**2 / (2 * sigma**2) ),
                     label = "std distribution",
                     linewidth=2, 
                     color='r' )

            plt.legend()
            filename = buildPlotFilename( options, key )
            plt.savefig( filename )

        E.info( "plotting P-value distribution" )
        
        key = "pvalue"
        plt.figure()

        x,bins,y = plt.hist( [r.pvalue for r in gat.iterator_results(annotator_results) ],
                             bins = numpy.arange( 0, 1.05, 0.025) ,
                             label = "pvalue" )

        plt.hist( [r.qvalue for r in gat.iterator_results(annotator_results) ],
                  bins = numpy.arange( 0, 1.05, 0.025) ,
                  label = "qvalue",
                  alpha=0.5 )

        plt.legend()

        # hist, bins = numpy.histogram( \
        #     [r.pvalue for r in gat.iterator_results(annotator_results) ],
        #     bins = 20 )
        # plt.plot( bins[:-1], hist, label = key )

        filename = buildPlotFilename( options, key )
        plt.savefig( filename )


    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
