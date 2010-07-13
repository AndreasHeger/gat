from cgat import *

import os, sys, re, optparse, collections, types
import numpy

import Bed
import IOTools
import Experiment as E

def readFromBed( filenames, name="track" ):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *name* and *contig*.
    *name* can "track" or "name".
    '''

    segment_lists = collections.defaultdict( lambda: collections.defaultdict(SegmentList))

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

class IntervalCollection(object):
    '''a collection of intervals.'''

    def __init__(self, name ):
        self.intervals = collections.defaultdict( lambda: collections.defaultdict(SegmentList))
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
        merged = collections.defaultdict( SegmentList )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                merged[contig].extend( segmentlist )
        self.intervals["merged"] = merged
        self.normalize()

    def sum( self ):
        '''remove all intervals not overlapping with intervals in other.'''
        s = 0
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                s += segmentlist.sum()
        return s

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
                    newlist = SegmentList( clone = segmentlist )
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

    def iteritems(self):
        return self.intervals.iteritems()
    
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
            for isochore in workspace.keys():
                val = counter( segs[isochore], annos[isochore], workspace[isochore] )
            f(track, annotation, aggregator(vals) )

    return counts

def iterator_results( annotator_results ):
    '''iterate over all results.'''
    for k1, v1 in annotator_results.iteritems():
        for k2, v2 in v1.iteritems():
            yield v2

def run( segments, 
         annotations, 
         workspace, 
         sampler, 
         counter,
         num_samples ):
    
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
        for x in xrange( num_samples ):
            E.info( "sampling: %s: %i/%i" % (track, x+1, num_samples))
            for isochore in segs.keys():
                # skip empty isochores
                if workspace[isochore].isEmpty: continue

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
    computeFDR( list(iterator_results(annotator_results)) )

    return annotator_results


