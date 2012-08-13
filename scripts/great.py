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

GREAT_RESULT = collections.defaultdict( "GREAT",
                                        ("segment",
                                         "annotation"
                                         "isochore",
                                         "nsegments_in_ws",
                                         "nsegments_in_anno",
                                         "nannotations_in_ws",
                                         "nannotations_in_segs",
                                         "coverage_annotations",
                                         "coverage_workspace",
                                         "pvalue" ) )

def main( argv ):

    tstart = time.time()

    segments, annotations, workspace = IO.buildSegments( options )

    E.info( "intervals loaded in %i seconds" % (time.time() - tstart) )

    # filter segments by workspace
    workspace = applyIsochores( segments, annotations, workspaces )

    ############################################
    ############################################
    ## per contig

    # compute percentage of base covered by annotations in workspace
    workspace2coverage = {}
    for contig, ww in workspace.iteritems():
        workspace2coverage[contig] = www.sum()

    annotation2bases = collections.defaultdict( dict )
    for annotation, aa in annotations.iteritems():
        for isochore, a in aa.iteritems():
            # need to truncate to workspace?
            annotation2coverage[annotation][isochore] = a.sum()
            
    results_per_contig = []

    # results per isochore 
    for segment, ss in segments.iteritems():
        for isochore, s in ss.iteritems():

            results = []
            temp = SegmentList( clone = ss )
            temp.intersect( workspace[isochore] )
            # number of segments in workspace
            nsegments_in_workspace = len(temp)
            
            for annotation, aa in annotations.iteritems():
                temp = SegmentList( clone = ss )
                temp.intersect( annotations[annotation][isochore] )
                # number of segments in annotation
                nsegments_in_annotation = len( temp )
                                         
                temp = SegmentList( clone = aa )
                temp.intersect( ss )
                nannotations_in_segments = len( temp )

                nannotations_in_workspace = len( aa )

                # p_A: proportion of bases covered by annotation
                coverage_annotation = annotation2coverage[annotation][isochore]
                coverage_workspace = workspace2coverage[isochore]
                percent_coverage_annotations = coverage_annotation / float( coverage_workspace )
                
                # GREAT binomial probability
                p_great = scipy.binom( ... )
                
                # hypergeometric probability - (annotations need to be non-overlapping!!!)
                p_hyper = hypergeom( nannotations_in_segments, 
                                     nannotations_in_workspace,
                                     nsegments_in_workspace )
                                    
                
                results.append( GREAT_RESULT._make( (
                            segment, annotation, isochore,
                            nsegments_in_workspace,
                            nsegments_in_annotations,
                            nannotations_in_workspace,
                            nannotations_in_segments,
                            coverage_annotations,
                            coverage_workspace,
                            percent_coverage_annotations,
                            p_great ) ) )
                
        results_per_contig.extend( results )


    # compute sums
    for segment in segments.tracks:
        for annotation in annotations.tracks:
            r = GREAT_RESULT._make( [ 
                    [ x.y for x in results_per_contig if x.segment == segment and x.annotation == annotation ]
                    for y in GREAT_RESULT.fields ] )
        

                            

