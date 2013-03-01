import unittest
import random, tempfile, shutil, os, re, gzip, sys, subprocess, collections
import numpy, math

from GatSegmentList import SegmentList
import GatEngine
import gat
import gat.IO

class Options: pass

class TestRunning( unittest.TestCase ):
    '''test running gat in various modes.
    '''

    sample_size = 10

    def setUp( self ):

        parser = gat.buildParser()

        options, args = parser.parse_args([])

        options.segment_files = ['data/segments.bed.gz']
        options.annotation_files = ['data/annotations.bed.gz']
        options.workspace_files = ['data/contigs.bed.gz' ]

        self.segments, self.annotations, workspaces, isochores = gat.IO.buildSegments( options )
        self.workspace = gat.IO.applyIsochores( self.segments, 
                                                self.annotations,
                                                workspaces,
                                                options,
                                                isochores )
        
        self.sampler = GatEngine.SamplerAnnotator( bucket_size = 1, nbuckets = 100000 )

        self.counters = [GatEngine.CounterNucleotideOverlap()]
        self.workspace_generator = GatEngine.UnconditionalWorkspace()
        
    def testNoMultiprocessing( self ):

        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size )
        
        self.assertEqual( len(annotator_results), 
                          len(self.segments) * len(self.annotations)
                          )

    def testOneCPU( self ):

        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size,
                                     num_threads = 1)
        
        self.assertEqual( len(annotator_results), 
                          len(self.segments) * len(self.annotations)
                          )


    def testTwoCPU( self ):
        return
        annotator_results = gat.run( self.segments,
                                     self.annotations,
                                     self.workspace,
                                     self.sampler,
                                     self.counters,
                                     workspace_generator = self.workspace_generator,
                                     num_samples = self.sample_size,
                                     num_threads = 2)
        
        self.assertEqual( len(annotator_results), 
                          len(self.segments) * len(self.annotations)
                          )

if __name__ == '__main__':
    unittest.main()
