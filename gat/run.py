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

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python script_template.py --help

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os, sys, re, optparse, collections, types
import numpy

import Experiment as E
import Bed
import IOTools
import gat

def readFromBed( filenames, name="track" ):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *name* and *contig*.
    *name* can "track" or "name".
    '''

    segment_lists = collections.defaultdict( lambda: collections.defaultdict(gat.SegmentList))

    if name == "track": f = lambda x: x.mTrack["name"]
    elif name == "name": f = lambda x: x.mFields[0]
    else: raise ValueError("unknown name: '%s'" %name )

    for filename in filenames:
        infile = IOTools.openFile( filename, "r")
        for bed in Bed.iterator( infile ):
            try:
                name = f(bed)
            except TypeError:
                name = "default"
            segment_lists[name][bed.contig].add( bed.start, bed.end )

    return segment_lists


class AnnotatorResult(object):
    
    format_observed = "%i"
    format_expected = "%5.2f"
    format_fold = "%6.4f"
    format_pvalue = "%6.4e"
    headers = ["track", "annotation", 
                "observed", "expected",
                "CI95low", "CI95high", 
                "fold", "pvalue" ]

    def __init__( self, 
                  track, annotation,
                  observed, samples ):
        self.track = track
        self.annotation = annotation
        self.observed = observed
        self.samples = numpy.array( samples, dtype= numpy.float )
        self.samples.sort()
        self.expected = numpy.mean(samples)
        self.fold = self.observed / self.expected
        self.stddev = numpy.std(samples)

        offset = int(0.05 * len(self.samples))
        self.lower95 = self.samples[ offset ]
        self.upper95 = self.samples[ -offset ]

        idx = numpy.searchsorted( self.samples, self.observed )
        self.pvalue = float(idx) / len(samples)
        # invert pvalue for test of depletion
        if self.fold > 1.0: self.pvalue = 1 - self.pvalue
        self.pvalue = max( self.pvalue, 1.0 / len(self.samples) )

    def __str__(self):

        if len(self.samples) < 10**6:
            format_pvalue = "%7.6f"
        else:
            format_pvalue = "%7.6e"

        return "\t".join( (self.track,
                           self.annotation,
                           self.format_observed % self.observed,
                           self.format_expected % self.expected,
                           self.format_expected % self.lower95,
                           self.format_expected % self.upper95,
                           self.format_expected % self.stddev,
                           self.format_fold % self.fold,
                           format_pvalue % self.pvalue ) )
        
class IntervalCollection(object):
    '''a collection of intervals.'''

    def __init__(self, name ):
        self.intervals = collections.defaultdict( lambda: collections.defaultdict(gat.SegmentList))
        self.name = name

    def load( self, filenames, name = "track" ):
        '''load segments from filenames and pre-process them.'''

        E.info( "%s: reading intervals from %i files" % (self.name, len(filenames)))
        self.intervals = readFromBed( filenames, name = name )
        E.info( "%s: read intervals for %i tracks" % (self.name, len(self.intervals) ))

    def normalize( self ):
        '''normalize segment lists individually.

        Remove empty contigs.
        '''

        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                if len(segmentlist) == 0: 
                    del vv[contig]
                else:
                    segmentlist.normalize()

    def outputStats( self, outfile ):
        '''output segment statistics.'''

        outfile.write( "section\ttrack\tcontig\tnsegments\tlength\n" )

        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % len(segmentlist), 
                         "%i" % segmentlist.sum() ) ) + "\n" )

    def merge( self ):
        '''merge all tracks into a single segment list
        creating a new track called 'merged`
        '''
        merged = collections.defaultdict( gat.SegmentList )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                merged[contig].extend( segmentlist )
        self.intervals["merged"] = merged
        self.normalize()

    def prune( self, other ):
        '''remove all intervals not overlapping with intervals in other.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.intersect( other[contig] )

    def sort( self ):
        '''sort all intervals lists.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.sort()

    def check( self ):
        '''check all intervals lists.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.sort()

    def restrict( self, restrict ):
        '''remove all tracks except those in restrict.'''
        if restrict in (list, tuple, set):
            r = set(restrict)
        else:
            r = set( [restrict,] )

        for track in self.intervals.keys():
            if track not in r:
                del self.intervals[track]

    def toIsochores( self, isochores ):
        '''split contigs in each track into isochores.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.items():
                for other_track, other_vv in isochores.iteritems():
                    newlist = gat.SegmentList( clone = segmentlist )
                    newlist.intersect( other_vv[contig] )
                    isochore = "%s.%s" % (contig, other_track)
                    vv[isochore] = newlist
                del vv[contig]

    @property
    def tracks(self): return self.intervals.keys()

    def keys(self): return self.intervals.keys()

    def add( self, track, contig, segmentlist ):
        self.intervals[track][contig] = segmentlist

    def __getitem__(self, key ):
        return self.intervals[key]
    
    def outputOverlapStats( self, outfile, other ):

        outfile.write( "section\ttrack\tcontig\toverlap\tlength\tdensity\n" )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                length = segmentlist.sum()
                if length == 0: continue
                overlap = segmentlist.overlapWithSegments(other[contig])
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % overlap,
                         "%i" % length,
                         "%f" % (float(overlap) / length)) ) + "\n" )

def computeCounts( counter, aggregator, 
                   segments, annotations, workspace,
                   append = False ):
    '''collect counts from *counter* between all combinations of *segments* and *annotations*.

    *aggregator* determines how values are combined between workspaces.

    If *append* is set, values for each track in *segmentns* are appended to a list for
    each track in *annotations*. Otherwise, a nested dictionary is returned.
    '''

    if append:
        counts = collections.defaultdict( list )
        f = lambda x,y,z: counts[y].append(z)
    else:
        counts = collections.defaultdict( lambda: collections.defaultdict( float ) )
        f = lambda x,y,z: counts[x].__setitem__( y, z )

    # collect counts per isochore
    for track in segments.tracks:
        segs = segments[track]
        for annotation in annotations.tracks:
            annos = annotations[annotation]
            vals = [ counter( segs[isochore], annos[isochore], workspace[isochore] ) for isochore in workspace.keys() ]
            f(track, annotation, aggregator(vals) )

    return counts
        
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
                      choices=("nucleotide-overlap", "nucleotide-density",),
                      help="quantity to test [default=%default]."  )

    parser.add_option("-n", "--num-samples", dest="num_samples", type="int", 
                      help="number of samples to compute [default=%default]."  )


    parser.set_defaults(
        annotation_files = [],
        segment_files = [],
        workspace_files = [],
        num_samples = 1000,
        nbuckets = 100000,
        bucket_size = 1,
        counter = "nucleotide-overlap",
        output_stats = "all",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ## TODO:
    ## * permit incremental runs
    ## * caching of samples, 
    ## * caching of segment size distribution
    ## * pre-compute workspace/segment intersection
    ## * parallel computation of samples 
    ## * speed

    ##################################################
    ##################################################
    ##################################################
    # process segments
    def dumpStats( coll, section ):
        if section in options.output_stats or "all" in options.output_stats:
            coll.outputStats( E.openOutputFile( section ) )
        
    # read one or more segment files
    segments = IntervalCollection( name = "segments " )
    segments.load( options.segment_files )
    dumpStats( segments, "stats_segments_raw" )
    segments.normalize()
    dumpStats( segments, "stats_segments_raw" )

    # read one or more annotations
    annotations = IntervalCollection( name = "annotations " )
    annotations.load( options.annotation_files )
    dumpStats( annotations, "stats_annotations_raw" )
    annotations.normalize()
    dumpStats( annotations, "stats_annotations_normed" )

    # read one or more workspaces
    workspaces = IntervalCollection( name = "workspaces " )
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
        isochores = IntervalCollection( name = "isochores" )
        isochores.load( options.isochore_files, name = "name" )
        dumpStats( isochores, "stats_isochores_raw" )

        # merge isochores and check if consistent (fully normalized)
        isochores.merge()
        isochores.sort()
        # discard individual isochores
        isochores.restrict( "merged" )
        # check that there are no overlapping segments within isochores
        isochores.check()

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
        counter = gat.CounterNucleotides()
    elif options.counter == "nucleotide-density":
        counter = gat.CounterDensity()
    else:
        raise ValueError("unknown counter '%s'" % options.counter )

    ##################################################
    ##################################################
    ##################################################
    # sample 
    # note: samples might need to be discarded if memory
    # is scarce.
    ##################################################
    samples = {}

    for track in segments.tracks:
        sample = IntervalCollection( name = track )
        segs = segments[track]
        for x in xrange( options.num_samples ):
            E.info( "sampling: %s: %i/%i" % (track, x+1, options.num_samples))
            for isochore in segs.keys():
                if workspace[isochore].isEmpty: 
                    E.warn( "skipping empty isochore: '%s' " % isochore)
                    continue

                r = sampler.sample( segs[isochore], workspace[isochore] )
                sample.add( "sample_%06i" % x, isochore, r )

        samples[track] = sample 

    ##################################################
    ##################################################
    ##################################################
    # collect counts from segments and for samples
    observed_counts = computeCounts( counter = counter,
                                     aggregator = sum,
                                     segments = segments,
                                     annotations = annotations,
                                     workspace = workspace )

    sampled_counts = {}
    for track in segments.tracks:
        sampled_counts[track] = computeCounts( counter = counter,
                                               aggregator = sum,
                                               segments = samples[track],
                                               annotations = annotations,
                                               workspace = workspace,
                                               append = True )
            
    ##################################################
    ##################################################
    ##################################################
    ## build annotator results
    ##################################################
    annotator_results = collections.defaultdict( dict )
    for track, r in observed_counts.iteritems():
        for annotation, observed in r.iteritems():
            annotator_results[track][annotation] = AnnotatorResult( \
                track = track,
                annotation = annotation,
                observed = observed,
                samples = sampled_counts[track][annotation] )

    ##################################################
    ##################################################
    ##################################################
    ## compute global fdr
    ##################################################
            

    ##################################################
    ##################################################
    ##################################################
    ## output
    ##################################################
    outfile = sys.stdout

    outfile.write("\t".join( AnnotatorResult.headers ) + "\n" )

    for track, r in annotator_results.iteritems():
        for annotation, result in r.iteritems():
            outfile.write( str(result) + "\n" )
            
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
