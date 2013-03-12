import unittest
import random, tempfile, shutil, os, re, gzip, sys, subprocess, collections
import numpy, math

from GatSegmentList import SegmentList
import GatEngine
import gat
import gat.IO
import gat.Stats

class Options: pass

class TestRunning( unittest.TestCase ):
    '''test running gat.

    Gat is being run in single-threaded and multi-threaded modes and the output
    is compared to reference gat results.

    The reference data has been created with the following commands::
       
        gat-run.py --segments=data/segments_multiple.bed.gz --annotations=data/annotations.bed.gz --workspace=data/workspace.bed.gz --num-samples=1000 > data/output_multiple.tsv
        gat-run.py --segments=data/segments_single.bed.gz --annotations=data/annotations.bed.gz --workspace=data/workspace.bed.gz --num-samples=1000 > data/output_single.tsv
    '''

    sample_size = 1000
    
    # similarity threshold for float comparisons
    threshold = 0.20

    filename_segments = [ 'data/segments_single.bed.gz' ]
    filename_annotations = [ 'data/annotations.bed.gz' ]
    filename_workspace = [ 'data/workspace.bed.gz' ]

    max_percent_difference = 10
    mean_percent_difference = 5

    def setUp( self ):

        parser = gat.buildParser()

        options, args = parser.parse_args([])

        options.segment_files = self.filename_segments
        options.annotation_files = self.filename_annotations
        options.workspace_files = self.filename_workspace

        self.segments, self.annotations, workspaces, isochores = gat.IO.buildSegments( options )
        self.workspace = gat.IO.applyIsochores( self.segments, 
                                                self.annotations,
                                                workspaces,
                                                options,
                                                isochores )
        
        self.sampler = GatEngine.SamplerAnnotator( bucket_size = 1, nbuckets = 100000 )

        self.counters = [GatEngine.CounterNucleotideOverlap()]
        self.workspace_generator = GatEngine.UnconditionalWorkspace()
        
        
        self.reference_data = gat.IO.readAnnotatorResults( 'data/output_single.tsv' )

    def compare( self, annotator_results ):

        self.assertEqual( len(annotator_results), 
                          len(self.segments) * len(self.annotations)
                          )

        self.assertEqual( len(annotator_results),
                          len(self.reference_data ) )

        self.assertEqual( sorted([x.track for x in annotator_results]),
                          sorted([x.track for x in self.reference_data]) )

        self.assertEqual( sorted([x.annotation for x in annotator_results]),
                          sorted([x.annotation for x in self.reference_data]) )
        
        # sort according to reference
        positions = dict( [ (y,x) for x,y in enumerate( [(x.track, x.annotation) for x in self.reference_data] )] )
        annotator_results.sort( key = lambda x: positions[(x.track,x.annotation)] )

        def _output_stats( s ):
            return ",".join( [ "%s=%s" % (x,y) for x,y in zip( summary.getHeaders(), 
                                                             str(summary).split("\t")) ] )
        for attr in ("expected", "fold", "pvalue" ):
            
            this = [ getattr(x,attr) for x in annotator_results ]
            ref = [ getattr(x,attr) for x in self.reference_data ]
        
            percent_differences = [ 100.0 * (x - y) / y for x,y in zip( this, ref ) ]
            summary = gat.Stats.Summary( percent_differences )
            mm = max(abs(summary.min), summary.max)
            self.assertLess( mm, self.max_percent_difference, 
                             "attribute '%s': maximum percent difference %f > %f: %s" \
                                 % (attr, 
                                    mm, self.max_percent_difference, 
                                 _output_stats( summary) ))

            self.assertLess( summary.mean, self.mean_percent_difference, 
                             "attribute '%s': mean percent difference %f > %f: %s" \
                                 % (attr, 
                                    summary.mean, self.max_percent_difference, 
                                 _output_stats( summary) ))
            

        # observed should be identical
        for ref in self.reference_data:
            this = [x for x in annotator_results if x.track == ref.track and x.annotation == ref.annotation ][0]
            self.assertEqual( ref.observed, this.observed )

    def testNoMultiprocessing( self ):

        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size )

        self.compare( annotator_results )

    def testOneCPU( self ):

        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size,
                                     num_threads = 1)
        
        self.compare( annotator_results )

    def testTwoCPU( self ):
        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size,
                                     num_threads = 2)
        
        self.compare( annotator_results )

if __name__ == '__main__':
    unittest.main()
