from cgat import *

import os, sys, re, optparse, collections, types, gzip, pprint
import numpy

import gat.Bed as Bed
import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats as Stats

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

class DummyAnnotatorResult:

    format_observed = "%i"
    format_expected = "%6.4f"
    format_fold = "%6.4f"
    format_pvalue = "%6.4e"

    def __init__( self ):
        pass

    @classmethod
    def _fromLine( cls, line ):
        x = cls()
        data = line[:-1].split("\t")
        x.track, x.annotation = data[:2]
        x.observed, x.expected, x.lower95, x.upper95, x.stddev, x.fold, x.pvalue, x.qvalue = \
            map(float, data[2:] )
        
        return x

    def __str__(self):
        return "\t".join( (self.track,
                           self.annotation,
                           self.format_observed % self.observed,
                           self.format_expected % self.expected,
                           self.format_expected % self.lower95,
                           self.format_expected % self.upper95,
                           self.format_expected % self.stddev,
                           self.format_fold % self.fold,
                           self.format_pvalue % self.pvalue,
                           self.format_pvalue % self.qvalue ) )


class UnconditionalSampler:

    def __init__(self,         
                 num_samples, samples, 
                 samples_outfile, 
                 outfile_sample_stats,
                 sampler,
                 workspace_generator, counter ):
        self.num_samples = num_samples
        self.samples = samples
        self.samples_outfile = samples_outfile
        self.sampler = sampler
        self.workspace_generator = workspace_generator
        self.counter = counter
        self.outfile_sample_stats = outfile_sample_stats
        
        if self.outfile_sample_stats:
            E.debug( "sample stats go to %s" % outfile_sample_stats)
            self.outfile_sample_stats.write( "sample\tisochore\tnsegments\tnnucleotides\tmean\tstd\tmin\tq1\tmedian\tq3\tmax\n" )

        self.last_sample_id = None
        self.all_lengths = []

    def outputSampleStats( self, sample_id, isochore, sample ):

        def _write( sample_id, isochore, lengths ):
            if self.outfile_sample_stats:
                self.outfile_sample_stats.write( "\t".join( map(str, (\
                                sample_id,
                                isochore,
                                len(lengths),
                                numpy.sum( lengths ),
                                numpy.mean( lengths ),
                                numpy.std( lengths ),
                                min( lengths ),
                                Stats.percentile( lengths, 0.25 ),
                                numpy.median( lengths ),
                                Stats.percentile( lengths, 0.75 ),
                                max(lengths) ) )) + "\n" )

            
        if sample_id != self.last_sample_id:
            if self.last_sample_id != None:
                _write( self.last_sample_id, "all", numpy.sort(numpy.array(self.all_lengths)) )
            self.all_lengths = []
            self.last_sample_id = sample_id

            
        if sample_id:
            l = sample.asLengths()
            self.all_lengths.extend( l  )
            _write( sample_id, isochore, numpy.sort( numpy.array( l )))
                        
    def sample( self, track, counts, counter, segs, annotations, workspace ):
        
        E.info( "performing unconditional sampling" )
                     
        counts_per_track = collections.defaultdict( list )

        E.info( "workspace without conditioning: %i segments, %i nucleotides" % \
                    (workspace.counts(),
                     workspace.sum() ) )

        temp_segs, _, temp_workspace = self.workspace_generator( segs, None, workspace )

        E.info( "workspace after conditioning: %i segments, %i nucleotides" % \
                     (workspace.counts(),
                      workspace.sum() ) )

        if workspace.sum() == 0:
            E.warn( "empty workspace - no computation performed" )
            return counts_per_track

        # compute samples unconditionally
        for x in xrange( self.num_samples ):
            # use textual sample ids to avoid parsing from dumped samples
            sample_id = str(x)
            E.debug( "progress: %s: %i/%i %i isochores" % (track, x+1, self.num_samples, len(segs.keys())))
            counts_per_isochore = collections.defaultdict( list )

            if self.samples_outfile: 
                self.samples_outfile.write("track name=%s\n" % sample_id)

            for isochore in segs.keys():
                counts.pairs += 1

                # skip empty isochores
                if temp_workspace[isochore].isEmpty or temp_segs[isochore].isEmpty: 
                    counts.skipped += 1
                    # too much information
                    # E.debug( "skipping empty isochore %s" % isochore )
                    continue

                # skip if read from cache
                if self.samples.hasSample( track, sample_id, isochore ): 
                    counts.loaded += 1
                    self.samples.load( track, sample_id, isochore )
                    r = self.samples[track][sample_id][isochore]
                    del self.samples[track][sample_id][isochore]
                else:
                    counts.sampled += 1
                    r = self.sampler.sample( temp_segs[isochore], temp_workspace[isochore] )
                    # too much information
                    # E.debug( "sample=%s, isochore=%s, segs=%i, sample=%i" % \
                    #             (sample_id, isochore, segs[isochore].sum(), r.sum()) )

                self.outputSampleStats( sample_id, isochore, r )

                # compute counts for each annotation/isochore and save
                for annotation in annotations.tracks:
                    annos = annotations[annotation]
                    counts_per_isochore[annotation].append( self.counter( r, 
                                                                          annos[isochore], 
                                                                          temp_workspace[isochore] ) )

                # save sample
                if self.samples_outfile: 
                    for start, end in r:
                        self.samples_outfile.write( "%s\t%i\t%i\n" % (isochore, start, end))

            # TODO: choose aggregator
            for annotation in annotations.tracks:
                counts_per_track[annotation].append( sum( counts_per_isochore[annotation] ) )
                
        self.outputSampleStats( None, "", [] )

        return counts_per_track

class ConditionalSampler( UnconditionalSampler ):
    
    def sample( self, track, counts, counter, segs, annotations, workspace ):
        '''conditional sampling - sample using only those 
        segments that contain both a segment and an annotation.

        return dictionary with counts per track
        '''

        E.info( "performing conditional sampling" )

        counts_per_annotation = collections.defaultdict( list )

        E.info( "workspace without conditioning: %i segments, %i nucleotides" % \
                     (workspace.counts(),
                      workspace.sum() ) )

        if workspace.sum() == 0:
            E.warn( "empty workspace - no computation performed" )
            return counts_per_track

        # compute samples conditionally
        for annoid, annotation in enumerate(annotations.tracks):

            annos = annotations[annotation]

            temp_segs, temp_annotations, temp_workspace = self.workspace_generator( segs, annos, workspace )

            E.info( "workspace for annotation %s: %i segments, %i nucleotides" % \
                        (annotation,
                         temp_workspace.counts(),
                         temp_workspace.sum() ) )

            for x in xrange( self.num_samples ):
                # use textual sample ids to avoid parsing from dumped samples
                sample_id = str(x)
                E.debug( "progress: %s: %i/%i %i/%i %i isochores" % (track, annoid+1, 
                                                                     len(annotations.tracks),
                                                                     x+1, self.num_samples, len(temp_segs.keys())))

                counts_per_isochore = []

                if self.samples_outfile: 
                    self.samples_outfile.write("track name=%s-%s\n" % (annotation, sample_id))

                for isochore in temp_segs.keys():
                    counts.pairs += 1

                    # skip empty isochores
                    if temp_workspace[isochore].isEmpty or temp_segs[isochore].isEmpty: 
                        counts.skipped += 1
                        continue

                    counts.sampled += 1
                    r = self.sampler.sample( temp_segs[isochore], temp_workspace[isochore] )
                    counts_per_isochore.append( self.counter( r, 
                                                              temp_annotations[isochore], 
                                                              temp_workspace[isochore] ) )


                    # save sample
                    if self.samples_outfile: 
                        for start, end in r:
                            self.samples_outfile.write( "%s\t%i\t%i\n" % (isochore, start, end))
                            
                    self.outputSampleStats( sample_id, isochore, r )

                # add sample to counts
                # TODO: choose aggregator
                sample_counts = sum( counts_per_isochore )
                counts_per_annotation[annotation].append( sample_counts )

            self.outputSampleStats( None, "", [] )

        return counts_per_annotation
    


###################################################################################
###################################################################################
###################################################################################
def sampleUnconditionally( track, counts, 
                           num_samples, samples, samples_outfile,
                           sampler,
                           workspace_generator,
                           counter, segs, annotations, workspace):
    '''unconditional sampling - sample over complete workspace.

    return dictionary with counts per annotation
    '''


    counts_per_track = collections.defaultdict( list )

    E.info( "workspace without conditioning: %i segments, %i nucleotides" % \
                (workspace.counts(),
                 workspace.sum() ) )

    temp_segs, _, temp_workspace = workspace_generator( segs, None, workspace )

    E.info( "workspace after conditioning: %i segments, %i nucleotides" % \
                 (workspace.counts(),
                  workspace.sum() ) )

    if workspace.sum() == 0:
        E.warn( "empty workspace - no computation performed" )
        return counts_per_track

    # compute samples unconditionally
    for x in xrange( num_samples ):
        # use textual sample ids to avoid parsing from dumped samples
        sample_id = str(x)
        E.debug( "progress: %s: %i/%i %i isochores" % (track, x+1, num_samples, len(segs.keys())))
        counts_per_isochore = collections.defaultdict( list )

        if samples_outfile: 
            samples_outfile.write("track name=%s\n" % sample_id)

        for isochore in segs.keys():
            counts.pairs += 1

            # skip empty isochores
            if temp_workspace[isochore].isEmpty or temp_segs[isochore].isEmpty: 
                counts.skipped += 1
                # too much information
                # E.debug( "skipping empty isochore %s" % isochore )
                continue

            # skip if read from cache
            if samples.hasSample( track, sample_id, isochore ): 
                counts.loaded += 1
                samples.load( track, sample_id, isochore )
                r = samples[track][sample_id][isochore]
                del samples[track][sample_id][isochore]
            else:
                counts.sampled += 1
                r = sampler.sample( temp_segs[isochore], temp_workspace[isochore] )
                # too much information
                # E.debug( "sample=%s, isochore=%s, segs=%i, sample=%i" % \
                #             (sample_id, isochore, segs[isochore].sum(), r.sum()) )

            # compute counts for each annotation/isochore and save
            for annotation in annotations.tracks:
                annos = annotations[annotation]
                counts_per_isochore[annotation].append( counter( r, annos[isochore], temp_workspace[isochore] ) )

            # save sample
            if samples_outfile: 
                for start, end in r:
                    samples_outfile.write( "%s\t%i\t%i\n" % (isochore, start, end))

        # TODO: choose aggregator
        for annotation in annotations.tracks:
            counts_per_track[annotation].append( sum( counts_per_isochore[annotation] ) )
            
    return counts_per_track



###################################################################################
###################################################################################
###################################################################################
def sampleConditionally( track, counts, 
                         num_samples, samples, samples_outfile,
                         sampler,
                         workspace_generator,
                         counter, segs, annotations, workspace):
    '''conditional sampling - sample using only those 
    segments that contain both a segment and an annotation.

    return dictionary with counts per track
    '''

    E.info( "performing conditional sampling" )

    counts_per_annotation = collections.defaultdict( list )

    E.info( "workspace without conditioning: %i segments, %i nucleotides" % \
                 (workspace.counts(),
                  workspace.sum() ) )
 
    # compute samples conditionally
    for annoid, annotation in enumerate(annotations.tracks):

        annos = annotations[annotation]

        temp_segs, temp_annotations, temp_workspace = workspace_generator( segs, annos, workspace )

        E.info( "workspace for annotation %s: %i segments, %i nucleotides" % \
                    (annotation,
                     temp_workspace.counts(),
                     temp_workspace.sum() ) )

        for x in xrange( num_samples ):
            # use textual sample ids to avoid parsing from dumped samples
            sample_id = str(x)
            E.debug( "progress: %s: %i/%i %i/%i %i isochores" % (track, annoid+1, 
                                                                 len(annotations.tracks),
                                                                 x+1, num_samples, len(temp_segs.keys())))

            counts_per_isochore = []

            if samples_outfile: 
                samples_outfile.write("track name=%s-%s\n" % (annotation, sample_id))

            for isochore in temp_segs.keys():
                counts.pairs += 1

                # skip empty isochores
                if temp_workspace[isochore].isEmpty or temp_segs[isochore].isEmpty: 
                    counts.skipped += 1
                    continue

                counts.sampled += 1
                r = sampler.sample( temp_segs[isochore], temp_workspace[isochore] )
                counts_per_isochore.append( counter( r, 
                                                     temp_annotations[isochore], 
                                                     temp_workspace[isochore] ) )


                # save sample
                if samples_outfile: 
                    for start, end in r:
                        samples_outfile.write( "%s\t%i\t%i\n" % (isochore, start, end))
            
            # add sample to counts
            # TODO: choose aggregator
            sample_counts = sum( counts_per_isochore )
            counts_per_annotation[annotation].append( sample_counts )

    return counts_per_annotation

def run( segments, 
         annotations, 
         workspace, 
         sampler, 
         counter,
         workspace_generator,
         **kwargs ):
    '''run an enrichment analysis.

    segments: an IntervalCollection
    workspace: an IntervalCollection
    annotations: an IntervalCollection

    kwargs recognized are:

    cache 
       filename of cache

    num_samples
       number of samples to compute

    output_counts
       output counts to filename 

    output_samples_pattern
       if given, output samles to these files, one per segment

    outfile_sample_stats
       output filename with sample stats

    sample_files
       if given, read samples from these files.

    fdr
       method to compute qvalues
    '''

    ## get arguments
    num_samples = kwargs.get( "num_samples", 10000 )
    cache = kwargs.get( "cache", None )
    output_counts = kwargs.get( "output_counts", None )
    sample_files = kwargs.get( "sample_files", [] )
    output_samples_pattern = kwargs.get( "output_samples_pattern", None )
    outfile_sample_stats = kwargs.get( "outfile_sample_stats", None )

    if output_samples_pattern != None: 
        if "%s" not in output_samples_pattern:
            raise ValueError( "output_samples_pattern should contain at least one '%s'")

    ##################################################
    ##################################################
    ##################################################
    # collect observed counts from segments
    E.info( "collecting observed counts" )
    observed_counts = computeCounts( counter = counter,
                                     aggregator = sum,
                                     segments = segments,
                                     annotations = annotations,
                                     workspace = workspace,
                                     workspace_generator = workspace_generator )
    E.info( "collecting observed densities" )

    ##################################################
    ##################################################
    ##################################################
    # sample and collect counts
    ##################################################
    if cache:
        E.info( "samples are cached in %s" % cache)
        samples = SamplesCached( filename = cache )
    elif sample_files:
        if not output_samples_pattern:
            raise ValueError( "require output_samples_pattern if loading samples from files" )
        # build regex
        regex = re.compile( re.sub("%s", "(\S+)", output_samples_pattern ) )
        E.info( "loading samples from %i files" % len(sample_files) )
        samples = SamplesFile( filenames = sample_files,
                               regex = regex ) 
    else:
        samples = Samples()

    sampled_counts = {}
    old_sampled_counts = {}
    
    counts = E.Counter()

    ntracks = len(segments.tracks)

    for ntrack, track in enumerate(segments.tracks):

        segs = segments[track]

        E.info( "sampling: %s: %i/%i" % (track, ntrack+1, ntracks))

        if output_samples_pattern and not sample_files:
            filename = re.sub("%s", track, output_samples_pattern )
            E.debug( "saving samples to %s" % filename)
            dirname = os.path.dirname( filename )
            if dirname and not os.path.exists( dirname ): os.makedirs( dirname )
            if filename.endswith(".gz"):
                samples_outfile = gzip.open( filename, "w" )                
            else:
                samples_outfile = open( filename, "w" )
        else:
            samples_outfile = None

        if workspace_generator.is_conditional:
            outer_sampler = ConditionalSampler( num_samples, 
                                          samples, 
                                          samples_outfile, 
                                          outfile_sample_stats,
                                          sampler,
                                          workspace_generator, 
                                          counter )
        else:
            outer_sampler = UnconditionalSampler( num_samples, 
                                            samples, 
                                            samples_outfile, 
                                            outfile_sample_stats,
                                            sampler,
                                            workspace_generator, 
                                            counter )


        counts_per_track = outer_sampler.sample( track, counts, counter, segs, annotations, workspace )

        if samples_outfile: samples_outfile.close()

        sampled_counts[track] = counts_per_track
        
        # old code, refactor into loop to save samples
        if 0:
            E.info( "sampling stats: %s" % str(counts))
            if track not in samples:
                E.warn( "no samples for track %s" % track )
                continue

            # clean up samples
            del samples[track]

    E.info( "sampling finished: %s" % str(counts) )

    ##################################################
    ##################################################
    ##################################################
    ## build annotator results
    ##################################################
    E.info( "computing PValue statistics" )
    annotator_results = collections.defaultdict( dict )
    for track, r in observed_counts.iteritems():
        for annotation, observed in r.iteritems():
            temp_segs, temp_annos, temp_workspace = workspace_generator( segments[track], 
                                                                         annotations[annotation], 
                                                                         workspace )

            annotator_results[track][annotation] = AnnotatorResult( \
                track = track,
                annotation = annotation,
                observed = observed,
                samples = sampled_counts[track][annotation],
                track_segments = temp_segs,
                annotation_segments = temp_annos,
                workspace = temp_workspace)

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

    return annotator_results


