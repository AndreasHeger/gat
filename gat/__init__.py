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

class SegmentsSummary:
    '''summarize segments in a workspace.
    '''
    header = ("all_segments",
              "all_nucleotides",
              "segments_overlapping_workspace",
              "nucleotides_overlapping_workspace",
              "segments_outside_workspace",
              "nucleotides_outside_workspace",
              "truncated_segments",
              "truncated_nucleotides",
              "density_workspace",
              "proportion_truncated_segments",
              "proportion_extending_nucleotides",
              "summary_all_segments",
              "summary_segments_overlapping_workspace",
              "summary_truncated_segments" )
                            

    def __init__(self):
        pass

    def update( self, segments, workspace ):
        '''compute summary statistics between two segments lists *segments* and *workspace*.
        This method computes:
            * number of segments/nucleotides
            * number of segments/nucleotides overlapping workspace
            * number of segments split by a boundary
            * number of nucleotides extending outside boundary
            * segment size distribution untruncated at boundaries (min, max, mean, media, q1, q3)
            * segment size distributino truncated at boundaries
        '''

        self.all_segments = len(segments)
        self.all_nucleotides = segments.sum()

        # build segments overlapping workspace
        segments_overlapping_workspace = csegmentlist.SegmentList( clone = segments )
        segments_overlapping_workspace.filter( workspace )

        # build segments truncated by workspace
        truncated_segments = csegmentlist.SegmentList( clone = segments_overlapping_workspace )
        truncated_segments.intersect( workspace )

        segments_extending_workspace = csegmentlist.SegmentList( clone = segments )
        segments_extending_workspace.subtract( segments_overlapping_workspace )

        # compute numbers
        self.segments_overlapping_workspace = len( truncated_segments )
        self.nucleotides_overlapping_workspace = truncated_segments.sum()

        self.segments_outside_workspace = self.all_segments - self.segments_overlapping_workspace
        self.nucleotides_outside_workspace = self.all_nucleotides - self.nucleotides_overlapping_workspace

        self.truncated_segments = len(segments_extending_workspace )
        self.truncated_nucleotides = segments_extending_workspace.sum()

        self.summary_all_segments = segments.summarize()
        self.summary_segments_overlapping_workspace = segments_overlapping_workspace.summarize()
        self.summary_truncated_segments = truncated_segments.summarize()

        workspace_size = workspace.sum()

        self.density_workspace, self.proportion_truncated_segments, self.proportion_extending_nucleotides = 0,0,0
        # some quality metrics
        if workspace_size > 0:
            self.density_workspace = float(self.nucleotides_overlapping_workspace) / workspace_size
        if self.segments_overlapping_workspace > 0:
            self.proportion_truncated_segments = float(self.truncated_segments) / self.segments_overlapping_workspace
            self.proportion_extending_nucleotides = float(self.truncated_nucleotides) / segments_overlapping_workspace.sum()

    def __str__(self):
        return "\t".join( [ "%i" % self.all_segments,
                            "%i" % self.all_nucleotides,
                            "%i" % self.segments_overlapping_workspace,
                            "%i" % self.nucleotides_overlapping_workspace,
                            "%i" % self.segments_outside_workspace,
                            "%i" % self.nucleotides_outside_workspace,
                            "%i" % self.truncated_segments,
                            "%i" % self.truncated_nucleotides,
                            "%f" % self.density_workspace,
                            "%f" % self.proportion_truncated_segments,
                            "%f" % self.proportion_extending_nucleotides,
                            "%s" % str(self.summary_all_segments),
                            "%s" % str(self.summary_segments_overlapping_workspace),
                            "%s" % str(self.summary_truncated_segments) ] )
                            

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
                        
    def sample( self, track, counts, counters, segs, annotations, workspace ):
        '''sample and return counts.

        Return a list of counted results for each counter.
        '''

        report_interval = 100

        # rebuild non-isochore annotations and workspace
        contig_annotations = annotations.clone()
        contig_annotations.fromIsochores()
        contig_annotations.name = "contig_" + annotations.name

        contig_workspace = workspace.clone()
        contig_workspace.fromIsochores()
        
        E.info( "performing unconditional sampling" )
        counts_per_track = [ collections.defaultdict( list ) for x in counters ]

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

        # work = [ (x, 
        #           track, 
        #           self.sampler,
        #           temp_segs, 
        #           annotations,
        #           contig_annotations,
        #           temp_workspace,
        #           contig_workspace,
        #           counts, counters,
        #           self.samples_outfile) for x in range(self.num_samples) ]

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

        n = len(work)

        E.info( "sampling started" )
        results = []
        if self.num_threads == 0:
            for i, w in enumerate(work):
                r = computeSample( w )
                if i % report_interval == 0:
                    E.info( "%i/%i done (%5.2f)" % (i, n, 100.0 * i / n ))
                results.append( r )
        else:
            E.info("setting up shared data structures")
            annotations.share()
            contig_annotations.share()
            contig_workspace.share( "contig_workspace" )
            temp_segs.share( "generated_segments" )
            temp_workspace.share( "generated_workspace" )
            
            E.info("generating processpool with %i threads for %i items" % (self.num_threads, len(work) ))
            pool = multiprocessing.Pool( self.num_threads )

            rs = pool.map_async(computeSample, work )

            for i, r in enumerate(pool.imap_unordered(computeSample, work)):
                if i % report_interval == 0:
                    E.info( "%i/%i done (%5.2f)" % (i, n, 100.0 * i / n ))
                results.append( r )

            pool.close()
            pool.join()

            # results = threadpool.map( computeSample, work )
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

        # This method needs to re-factored to remove isochores before counting
        # and to work with multiple counters
        # also: parallel computation
        raise NotImplementedError

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
    # computing summary statistics for segments and
    # annotations
    if "segment_metrics" in outfiles:
        outfile = outfiles["segment_metrics"]
        E.info( "computing summary metrics for segments" )
    
        outfile.write( "metric\t%s\n" % "\t".join(Stats.Summary().getHeaders() ))
        observed_counts = []
        for ntrack, track in enumerate(segments.tracks):
            segs = segments[track]
            stats_per_isochore = []
            for isochore, ss in segs.iteritems():
                stats = SegmentsSummary()
                stats.update( ss, workspace[isochore] )
                stats_per_isochore.append( stats )

            for attribute in ("all_segments", "all_nucleotides",
                              "segments_overlapping_workspace",
                              "nucleotides_overlapping_workspace",
                              "nucleotides_outside_workspace",
                              "truncated_segments",
                              "truncated_nucleotides",
                              "density_workspace",
                              "proportion_truncated_segments",
                              "proportion_extending_nucleotides",
                              ):
                values = [ getattr( x, attribute ) for x in stats_per_isochore ]
                outfile.write("%s\t%s\n" % (attribute, Stats.Summary( values )) )
        outfile.flush()
        E.info( "wrote summary metrics for segments to %s" % str(outfile))

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
