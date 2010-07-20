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

    output_counts
       output counts to filename 
    '''

    ## get arguments
    num_samples = kwargs.get( "num_samples", 10000 )
    cache = kwargs.get( "cache", None )
    output_counts = kwargs.get( "output_counts", None )
    
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

    ##################################################
    ##################################################
    ##################################################
    ## build annotator results
    ##################################################
    E.info( "computing PValue statistics" )
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
    ## dump large table with counts
    ##################################################
    if output_counts:
        E.info( "writing counts to %s" % output_counts )
        output = list( iterator_results( annotator_results ) )
        outfile = open( output_counts, "w")
        outfile.write("sampleid" )
        for o in output:
            outfile.write("\t%s-%s" % (o.track, o.annotation) )
        outfile.write("\n")

        outfile.write("observed\t%s\n" % "\t".join(map(str, [o.observed for o in output ] ) ) )

        for x in xrange(num_samples):
            outfile.write( "%i\t%s\n" % \
                               (x, 
                                "\t".join(map(str, [o.getSample(x) for o in output ] ) ) ) )

    ##################################################
    ##################################################
    ##################################################
    ## compute global fdr
    ##################################################
    E.info( "computing FDR statistics" )
    computeFDR( list(iterator_results(annotator_results)) )

    return annotator_results

def fromCounts( filename ):
    '''build annotator results from a tab-separated table
    with counts.'''

    annotator_results = collections.defaultdict( dict )

    with open( filename, "r") as infile:

        E.info( "loading data")

        headers = infile.readline()[:-1].split("\t")[1:]
        observed = numpy.array( infile.readline()[:-1].split("\t")[1:], dtype = numpy.float)
        samples = numpy.loadtxt( infile, dtype=numpy.float, delimiter="\t" )

        E.info( "computing PValue statistics" )

        for x,header in enumerate(headers):
            track, annotation = header.split("-")
            annotator_results[track][annotation] = AnnotatorResult( \
                track = track,
                annotation = annotation,
                observed = observed[x],
                samples = samples[:,x+1] )

    ##################################################
    ##################################################
    ##################################################
    ## compute global fdr
    ##################################################
    E.info( "computing FDR statistics" )
    # computeFDR( list(iterator_results(annotator_results)) )

    return annotator_results


