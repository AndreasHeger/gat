from cgat import *

import os, sys, re, optparse, collections, types
import numpy

import Bed
import IOTools
import Experiment as E

def readFromBedOld( filenames, name = "track" ):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *contig* and *track*.
    
    If no track is given, the *name* attribute is taken.
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
         **kwargs ):
    '''run an enrichment analysis.

    kwargs recognized are:

    cache 
       filename of cache

    num_samples
       number of samples to compute

    '''

    ## get arguments
    num_samples = kwargs.get( "num_samples", 10000 )
    cache = kwargs.get( "cache", None )

    ##################################################
    ##################################################
    ##################################################
    # collect observed counts from segments
    E.info( "collecting observed counts" )
    observed_counts = computeCounts( counter = counter,
                                     aggregator = sum,
                                     segments = segments,
                                     annotations = annotations,
                                     workspace = workspace )

    ##################################################
    ##################################################
    ##################################################
    # sample and collect counts
    ##################################################
    E.info( "starting sampling" )

    if cache:
        E.info( "samples are cached in %s" % cache)
        samples = SamplesCached( filename = cache )
    else:
        samples = Samples( )

    sampled_counts = {}
    
    for track in segments.tracks:
        segs = segments[track]
        E.info( "sampling: %s" % (track))
        for x in xrange( num_samples ):
            E.debug( "progress: %s: %i/%i" % (track, x+1, num_samples))
            for isochore in segs.keys():
                # skip empty isochores
                if workspace[isochore].isEmpty or segs[isochore].isEmpty: continue
                # skip if read from cache
                if samples.hasSample( track, x, isochore ): 
                    samples.load( track, x, isochore )
                else:
                    r = sampler.sample( segs[isochore], workspace[isochore] )
                    samples.add( track, x, isochore, r )

        sampled_counts[track] = computeCounts( counter = counter,
                                               aggregator = sum,
                                               segments = samples[track],
                                               annotations = annotations,
                                               workspace = workspace,
                                               append = True )
        # clean up samples
        del samples[track]

    E.info( "sampling finished" )

    E.info( "computing PValue statistics" )
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
    E.info( "computing FDR statistics" )
    computeFDR( list(iterator_results(annotator_results)) )

    return annotator_results


