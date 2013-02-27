from cgat import *

import os, sys, re, optparse, collections, types, gzip, pprint, time
import numpy

import gat.Bed as Bed
import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats as Stats

import multiprocessing.pool

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
        x.counter = "na"
        x.observed, x.expected, x.lower95, x.upper95, x.stddev, x.fold, x.pvalue, x.qvalue = \
            map(float, data[2:10] )
        if len(data) > 10:
            (track_nsegments,
            track_size,
            track_density,
            annotation_nsegments,
            annotation_size,
            annotation_density,
            overlap_nsegments,
            overlap_size,
            overlap_density,
            percent_overlap_nsegments_track,
            percent_overlap_size_track,
            percent_overlap_nsegments_annotation,
            percent_overlap_size_annotation) = map(float, data[10:23])
            

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


def computeSample( args ):
    '''perform a simple sample

    Unconditional sampling.
    '''
    ( sample_id, 
      track, 
      sampler, 
      segs, 
      annotations, 
      contig_annotations,
      workspace, 
      contig_workspace,
      counters,
      samples_outfile) = args

    counts = Experiment.Counter()

    sample_id = str(sample_id)
    counts_per_isochore = collections.defaultdict( list )

    if samples_outfile: 
        samples_outfile.write("track name=%s\n" % sample_id)

    sample = IntervalDictionary()

    for isochore in segs.keys():

        counts.pairs += 1

        # skip empty isochores
        if workspace[isochore].isEmpty or segs[isochore].isEmpty: 
            counts.skipped += 1
            continue

        counts.sampled += 1
        r = sampler.sample( segs[isochore], workspace[isochore] )

        # TODO : activate
        # self.outputSampleStats( sample_id, isochore, r )

        sample.add( isochore, r )

        # save sample
        if samples_outfile: 
            for start, end in r:
                samples_outfile.write( "%s\t%i\t%i\n" % (isochore, start, end))

    # re-combine isochores
    # adjacent intervals are merged.
    sample.fromIsochores()

    counts_per_track = [ collections.defaultdict(float) for x in counters ]
    # compute counts for each counter
    for counter_id, counter in enumerate(counters):
        # TODO: choose aggregator
        for annotation in annotations.tracks:
            counts_per_track[counter_id][annotation] = sum( [
                    counter( sample[contig],
                             contig_annotations[annotation][contig],
                             contig_workspace[contig])
                    for contig in sample.keys() ] )

    return counts_per_track

class UnconditionalSampler:

    def __init__(self,         
                 num_samples, samples, 
                 samples_outfile, 
                 outfile_sample_stats,
                 sampler,
                 workspace_generator, 
                 counters,
                 num_threads = 1 ):
        self.num_samples = num_samples
        self.samples = samples
        self.samples_outfile = samples_outfile
        self.sampler = sampler
        self.workspace_generator = workspace_generator
        self.counters = counters
        self.outfile_sample_stats = outfile_sample_stats
        
        if self.outfile_sample_stats:
            E.debug( "sample stats go to %s" % outfile_sample_stats)
            self.outfile_sample_stats.write( "sample\tisochore\tnsegments\tnnucleotides\tmean\tstd\tmin\tq1\tmedian\tq3\tmax\n" )

        self.last_sample_id = None
        self.all_lengths = []
        self.num_threads = num_threads

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

    def computeSamples( self, work, report_interval = 100 ):
        '''compute samples according to work.

        returns a list of results.
        '''
        
        n = len(work)

        results = []

        if self.num_threads == 0:
            for i, w in enumerate(work):
                r = computeSample( w )
                if i % report_interval == 0:
                    E.info( "%i/%i done (%5.2f)" % (i, n, 100.0 * i / n ))
                results.append( r )
        else:
            E.info("generating processpool with %i threads for %i items" % (self.num_threads, len(work) ))
            pool = multiprocessing.Pool( self.num_threads )

            rs = pool.map_async(computeSample, work )

            for i, r in enumerate(pool.imap_unordered(computeSample, work)):
                if i % report_interval == 0:
                    E.info( "%i/%i done (%5.2f)" % (i, n, 100.0 * i / n ))
                results.append( r )

            pool.close()
            pool.join()

        return results

    def sample( self, track, counts, counters, segs, annotations, workspace ):
        '''sample and return counts.

        Return a list of counted results for each counter.
        '''

        E.info( "performing unconditional sampling" )
        counts_per_track = [ collections.defaultdict( list ) for x in counters ]

        # rebuild non-isochore annotations and workspace
        contig_annotations = annotations.clone()
        contig_annotations.fromIsochores()
        contig_annotations.name = "contig_" + annotations.name

        contig_workspace = workspace.clone()
        contig_workspace.fromIsochores()

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

        work = [ (x, 
                  track, 
                  self.sampler,
                  temp_segs, 
                  annotations,
                  contig_annotations,
                  temp_workspace,
                  contig_workspace,
                  counters, 
                  self.samples_outfile
                  ) for x in range(self.num_samples) ]

    
        E.info("setting up shared data for multi-processing")
        annotations.share()
        contig_annotations.share()
        contig_workspace.share( "contig_workspace" )
        temp_segs.share( "generated_segments" )
        temp_workspace.share( "generated_workspace" )

        E.info( "sampling started" )
        results = self.computeSamples( work )
        E.info( "sampling completed" )

        # collate results
        for result in results:
            for counter_id, counter in enumerate(counters):
                for annotation in annotations.tracks:
                    counts_per_track[counter_id][annotation].append(
                        result[counter_id][annotation] )

        self.outputSampleStats( None, "", [] )

        return counts_per_track

class ConditionalSampler( UnconditionalSampler ):
    
    def sample( self, track, counts, counters, segs, annotations, workspace ):
        '''conditional sampling - sample using only those 
        segments that contain both a segment and an annotation.

        return dictionary with counts per track
        '''

        E.info( "performing conditional sampling" )
        counts_per_track = [ collections.defaultdict( list ) for x in counters ]

        # rebuild non-isochore annotations and workspace
        contig_annotations = annotations.clone()
        contig_annotations.fromIsochores()
        contig_annotations.name = "contig_" + annotations.name

        contig_workspace = workspace.clone()
        contig_workspace.fromIsochores()

        E.info("setting up shared data for multi-processing")
        annotations.share()
        contig_annotations.share()
        contig_workspace.share( "contig_workspace" )

        E.info( "workspace without conditioning: %i segments, %i nucleotides" % \
                     (workspace.counts(),
                      workspace.sum() ) )

        if workspace.sum() == 0:
            E.warn( "empty workspace - no computation performed" )
            return counts_per_track

        # compute samples conditionally - need to proceed by annotation
        for annoid, annotation in enumerate(annotations.tracks):

            annos = annotations[annotation]

            temp_segs, temp_annotations, temp_workspace = self.workspace_generator( segs, annos, workspace )
            
            # set up sharing
            temp_segs.share( "generated_segments" )
            temp_workspace.share( "generated_workspace" )

            E.info( "workspace for annotation %s: %i segments, %i nucleotides" % \
                        (annotation,
                         temp_workspace.counts(),
                         temp_workspace.sum() ) )

            work = [ (annoid, 
                      track, 
                      self.sampler,
                      temp_segs, 
                      annotations,
                      contig_annotations,
                      temp_workspace,
                      contig_workspace,
                      counters, 
                      self.samples_outfile
                      ) for x in range(self.num_samples) ]
            
            E.info( "sampling for annotation '%s' started" % annotation)
            results = self.computeSamples( work )
            E.info( "sampling for annotation '%s' completed" % annotation)

            for result in results:
                for counter_id, counter in enumerate(counters):
                    counts_per_track[counter_id][annotation].append( result[counter_id][annotation] )

        return counts_per_track

def run( segments, 
         annotations, 
         workspace, 
         sampler, 
         counters,
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

    output_counts_pattern
       output counts to filename 

    output_samples_pattern
       if given, output samles to these files, one per segment

    sample_files
       if given, read samples from these files.

    fdr
       method to compute qvalues

    outfiles
       dictionary of optional additional output files.

    pseudo_count
       pseudo_count to add to observed and expected values

    reference
       data with reference observed and expected values.
    '''

    ## get arguments
    num_samples = kwargs.get( "num_samples", 10000 )
    cache = kwargs.get( "cache", None )
    output_counts_pattern = kwargs.get( "output_counts_pattern", None )
    sample_files = kwargs.get( "sample_files", [] )
    pseudo_count = kwargs.get( "pseudo_count", 1.0 )
    reference = kwargs.get( "reference", None )
    output_samples_pattern = kwargs.get( "output_samples_pattern", None )
    outfiles = kwargs.get( "outfiles", {} )
    num_threads = kwargs.get( "num_threads", 1 )

    ##################################################
    ##################################################
    ##################################################
    # computing summary metrics for segments
    if "segment_metrics" in outfiles:
        IO.outputMetrics(  outfiles["segment_metrics"],
                           segments, workspace )

    ##################################################
    ##################################################
    ##################################################
    # collect observed counts from segments
    E.info( "collecting observed counts" )
    observed_counts = []
    for counter in counters:
        observed_counts.append( computeCounts( counter = counter,
                                               aggregator = sum,
                                               segments = segments,
                                               annotations = annotations,
                                               workspace = workspace,
                                               workspace_generator = workspace_generator ) )

    ##################################################
    ##################################################
    ##################################################
    # sample and collect counts
    ##################################################
    E.info( "starting sampling" )

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
                                                outfiles.get( "sample_stats", None ),
                                                sampler,
                                                workspace_generator, 
                                                counters,
                                                num_threads = num_threads )
        else:
            outer_sampler = UnconditionalSampler( num_samples, 
                                                  samples, 
                                                  samples_outfile, 
                                                  outfiles.get( "sample_stats", None ),
                                                  sampler,
                                                  workspace_generator, 
                                                  counters,
                                                  num_threads = num_threads )

        counts_per_track = outer_sampler.sample( track, counts, counters, segs, annotations, workspace )

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

    annotator_results = list()

    counter_id = 0
    for counter, observed_count in zip( counters, observed_counts ):

        for track, r in observed_count.iteritems():
            for annotation, observed in r.iteritems():
                temp_segs, temp_annos, temp_workspace = workspace_generator( segments[track], 
                                                                             annotations[annotation], 
                                                                             workspace )

                # if reference is given, p-value will indicate difference
                # The test that track and annotation are present is done elsewhere
                if reference:
                    ref = reference[track][annotation]
                else:
                    ref = None
                annotator_results.append( AnnotatorResultExtended( \
                        track = track,
                        annotation = annotation,
                        counter = counter.name,
                        observed = observed,
                        samples = sampled_counts[track][counter_id][annotation],
                        track_segments = temp_segs,
                        annotation_segments = temp_annos,
                        workspace = temp_workspace,
                        reference = ref,
                        pseudo_count = pseudo_count ) )
        counter_id += 1

    ##################################################
    ##################################################
    ##################################################
    ## dump large table with counts
    ##################################################
    if output_counts_pattern:
        for counter in counters:
            name = counter.name
            filename = re.sub("%s", name, output_counts_pattern )
            
            E.info( "writing counts to %s" % filename )
            output = [ x for x in annotator_results if x.counter == name ]
            outfile = IOTools.openFile( filename, "w")
            outfile.write("track\tannotation\tobserved\tcounts\n" )
        
            for o in output:
                outfile.write( "%s\t%s\t%i\t%s\n" % \
                                   (o.track, o.annotation,
                                    o.observed,
                                    ",".join(["%i" % x for x in o.samples]) ) )
                
    return annotator_results

def fromCounts( filename ):
    '''build annotator results from a tab-separated table
    with counts.'''

    annotator_results = []

    with IOTools.openFile( filename, "r") as infile:

        E.info( "loading data")

        header = infile.readline()
        if not header== "track\tannotation\tobserved\tcounts\n":
            raise ValueError("%s not a counts file: got %s" % (infile, header) )
        
        for line in infile:
            track, annotation, observed, counts = line[:-1].split( "\t" )
            samples = numpy.array( map(float, counts.split(",")), dtype=numpy.float )
            observed = float(observed)
            annotator_results.append( gat.AnnotatorResult( 
                track = track,
                annotation = annotation,
                counter = "na",
                observed = observed,
                samples = samples ) )

    return annotator_results
