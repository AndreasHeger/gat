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

import os, sys, re, optparse, collections, types
import numpy

import Experiment as E

import IOTools
import gat

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--annotation-file", dest="annotation_files", type="string", action="append",
                      help="filename with annotations [default=%default]."  )

    parser.add_option("-s", "--segment-file", dest="segment_files", type="string", action="append",
                      help="filename with segments [default=%default]."  )

    parser.add_option("-w", "--workspace-file", dest="workspace_files", type="string", action="append",
                      help="filename with workspace segments [default=%default]."  )

    parser.add_option("-i", "--isochore-file", dest="isochore_files", type="string", action="append",
                      help="filename with isochore segments [default=%default]."  )

    parser.add_option("-c", "--counter", dest="counter", type="choice",
                      choices=("nucleotide-overlap", 
                               "nucleotide-density",
                               "segment-overlap", ),
                      help="quantity to test [default=%default]."  )

    parser.add_option("-n", "--num-samples", dest="num_samples", type="int", 
                      help="number of samples to compute [default=%default]."  )

    parser.add_option("-o", "--order", dest="output_order", type="choice",
                      choices = ( "track", "annotation", "fold", "pvalue", "qvalue" ),
                      help="output order of results [default=%default]."  )


    parser.set_defaults(
        annotation_files = [],
        segment_files = [],
        workspace_files = [],
        num_samples = 1000,
        nbuckets = 100000,
        bucket_size = 1,
        counter = "nucleotide-overlap",
        output_stats = "all",
        output_order = "fold",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

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
    # process segments
    def dumpStats( coll, section ):
        if section in options.output_stats or "all" in options.output_stats:
            coll.outputStats( E.openOutputFile( section ) )
        
    # read one or more segment files
    segments = gat.IntervalCollection( name = "segments " )
    segments.load( options.segment_files )
    dumpStats( segments, "stats_segments_raw" )
    segments.normalize()
    dumpStats( segments, "stats_segments_raw" )

    # read one or more annotations
    annotations = gat.IntervalCollection( name = "annotations " )
    annotations.load( options.annotation_files )
    dumpStats( annotations, "stats_annotations_raw" )
    annotations.normalize()
    dumpStats( annotations, "stats_annotations_normed" )

    # read one or more workspaces
    workspaces = gat.IntervalCollection( name = "workspaces " )
    workspaces.load( options.workspace_files )
    dumpStats( workspaces, "stats_workspaces_raw" )
    workspaces.normalize()
    dumpStats( workspaces, "stats_workspaces_normed" )

    # intersect workspaces to build a single workspace
    workspaces.merge()
    dumpStats( workspaces, "stats_workspaces_merged" )

    # use merged workspace only, discard others
    workspaces.restrict("merged")

    if options.isochore_files:

        # read one or more isochore files
        isochores = gat.IntervalCollection( name = "isochores" )
        isochores.load( options.isochore_files, name = "name" )
        dumpStats( isochores, "stats_isochores_raw" )

        # merge isochores and check if consistent (fully normalized)
        isochores.sort()

        # check that there are no overlapping segments within isochores
        isochores.check()
        # TODO: flag is_normalized not properly set
        isochores.normalize()

        # check that there are no overlapping segments between isochores

        # intersect isochores and workspaces, segments and annotations
        workspaces.toIsochores( isochores )
        annotations.toIsochores( isochores )
        segments.toIsochores( isochores )

        dumpStats( workspaces, "stats_workspaces_isochores" )
        dumpStats( annotations, "stats_annotations_isochores" )
        dumpStats( segments, "stats_segments_isochores" )
    
    workspace = workspaces["merged"] 

    # prune segments and annotations keeping only
    # those overlapping the workspace
    #segments.prune( workspace )
    #segments.outputStats( options.stdout )
    #annotations.prune( workspace )
    #annotations.outputStats( options.stdout )

    # output segment densities per workspace
    for track in segments.tracks:
        workspaces.outputOverlapStats( options.stdout, segments[track] )

    ##################################################
    ##################################################
    ##################################################
    # initialize sampler
    sampler = gat.SamplerAnnotator(
        bucket_size = options.bucket_size,
        nbuckets = options.nbuckets )

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
    ## compute
    ##################################################
    annotator_results = gat.run( segments, 
                                 annotations, 
                                 workspace,
                                 sampler, 
                                 counter,
                                 num_samples = options.num_samples )
    
    ##################################################
    ##################################################
    ##################################################
    ## output
    ##################################################
    outfile = sys.stdout

    outfile.write("\t".join( gat.AnnotatorResult.headers ) + "\n" )

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
        outfile.write( str(result) + "\n" )
            
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
