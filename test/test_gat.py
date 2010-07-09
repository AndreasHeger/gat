"""test high-level interface."""

import unittest
import random, tempfile, shutil, os, re, gzip
import gat
import numpy

class TestSegmentList( unittest.TestCase ):
    
    def testCreateAndClear( self ):
        s = gat.SegmentList()
        self.assertEqual( 0, len(s) )
        s.add( 0, 100)
        self.assertEqual( 1, len(s) )
        s.clear()
        self.assertEqual( 0, len(s) )

    def testNormalize1( self ):
        '''non-overlapping segments.'''

        ss = [ (x, x + 10 ) for x in range( 0, 1000, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()

        self.assertEqual( len(s), 10 )
        self.assertEqual( s.sum(), 100 )

    def testNormalize2( self ):
        '''overlapping segments.'''

        ss = [ (x, x + 1000 ) for x in range( 0, 1000, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()

        self.assertEqual( len(s), 1 )
        self.assertEqual( s.sum(), 1000 )

    def testNormalize3( self ):
        '''non-overlapping but adjacent segments.'''

        ss = [ (x, x + 100 ) for x in range( 0, 1000, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()

        self.assertEqual( len(s), 10 )
        self.assertEqual( s.sum(), 1000 )

    def testNormalize4( self ):
        # test multiple interleaved segments
        ss = [ (x, x + 100 ) for x in range( 0, 1000, 10) ]
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()
        self.assertEqual( len(s), 1 )
        self.assertEqual( s.sum(), 1090 )

    def testExtend( self ):
        
        s1 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 0, 1000, 100) ] )
        s2 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 2000, 3000, 100) ] )
        s1.extend(s2 )
        self.assertEqual( s1.sum(), s2.sum() * 2 )
        self.assertEqual( len(s1), len(s2) * 2 )

class TestSegmentListOverlap( unittest.TestCase ):
    
    def setUp( self ):
        self.a = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 0, 1000, 100) ), normalize = True )

    def testOverlapFull( self ):
        self.assertEqual( self.a.overlapWithRange( 0, 1000), self.a.sum() )

    def testOverlapHalf( self ):
        self.assertEqual( self.a.overlapWithRange( 0, 500), self.a.sum() / 2)
        self.assertEqual( self.a.overlapWithRange( 500, 1000), self.a.sum() / 2)

    def testOverlapAfter( self ):
        self.assertEqual( self.a.overlapWithRange( 900, 910), 10 )
        self.assertEqual( self.a.overlapWithRange( 905, 915), 5 )

    def testOverlapNone( self ):
        self.assertEqual( self.a.overlapWithRange( 1000, 2000), 0 )
        self.assertEqual( self.a.overlapWithRange( 2000, 3000), 0 )

    def testOverlapAll( self ):
        for x in range( 0, 1000, 100):
            for y in range( 0, 10 ):
                self.assertEqual( self.a.overlapWithRange( x+y,x+y+1), 1, \
                                      "no overlap failure at %i: %i" % (x+y, self.a.overlapWithRange( x+y,x+y+1)))
            for y in range( 10, 100 ):
                self.assertEqual( self.a.overlapWithRange( x+y,x+y+1), 0, "overlap failure at %i" % (x+y) )

class TestSegmentListIntersection( unittest.TestCase):

    def setUp( self ):
        self.a = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 0, 1000, 100) ), normalize = True )

    def testIntersectionFull( self ):
        b = gat.SegmentList( iter = [ (0, 1000) ], normalize = True  ) 
        r = b.intersect( self.a )
        self.assertEqual( r.asList(), self.a.asList() )

    def testIntersectionSelf( self ):
        r = self.a.intersect( self.a )
        self.assertEqual( r.asList(), self.a.asList() )

    def testIntersectionCopy( self ):
        b = gat.SegmentList( clone = self.a )
        r = b.intersect( self.a )
        self.assertEqual( r.asList(), self.a.asList() )
        
    def testNoIntersection( self ):
        b = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 10, 1000, 100) ), normalize = True )
        r = b.intersect( self.a )
        self.assertEqual( r.asList(), [] )
        self.assertEqual( r.isEmpty(), True )

    def testPartialIntersection( self ):
        b = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 5, 1000, 100) ), normalize = True )
        r = b.intersect( self.a )
        self.assertEqual( len(r), len(self.a) )
        self.assertEqual( r.sum(), self.a.sum() / 2 )

class TestSamplerLength( unittest.TestCase ):

    ntests = 1000
    nsegments = 10000
    sampler = None

    def setUp( self ):
        
        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        self.segments = gat.SegmentList( iter = [ (x, x + numpy.random.randn() * 10.0 + 100.0  ) \
                                                      for x in range(0, 1000 * self.nsegments, 1000) ],
                                         normalize = True )

        self.histogram = self.segments.getLengthDistribution( 1, 1000 * self.nsegments )

    def testSampling( self ):

        if not self.sampler: return

        hs = self.sampler( self.histogram, 1 )

        samples = [hs.sample() for x in range(0, 1 * self.nsegments )]

        self.assertAlmostEqual( numpy.mean(samples),
                                100.0,
                                places = 0 )

        self.assertAlmostEqual( numpy.std(samples),
                                10.0,
                                places = 0 )
    
class TestSamplerLengthFast( TestSamplerLength ):

    sampler = gat.HistogramSampler

class TestSamplerLengthSlow( TestSamplerLength ):

    sampler = gat.HistogramSamplerSlow

class TestSamplerAnnotator( unittest.TestCase ):

    ntests = 1000

    def testTestSamplingSimple( self ):
        '''test if we sample the exactly right amount of nucleotides.'''
        nsegments = 10

        workspace = gat.SegmentList( iter = ( (x, x + 1000 ) \
                                                  for x in range( 0, 1000 * nsegments, 1000) ), 
                                     normalize = True )

        segments = gat.SegmentList( iter = ( (x,  x + 100)   \
                                        for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )
        
        sampler = gat.SamplerAnnotator()
        
        counts = numpy.zeros( 1000 * nsegments, numpy.int )

        for x in range( self.ntests):
            sample = sampler.sample( segments, workspace )
            self.assertEqual( sample.sum(), segments.sum() )

    def testTestSampling( self ):
        '''test if we sample the exactly right amount of nucleotides
        and the density is as expected.'''
        nsegments = 10000

        workspace = gat.SegmentList( iter = ( (x, x + 1000 ) \
                                                  for x in range( 0, 1000 * nsegments, 1000) ), 
                                     normalize = True )
        segments = gat.SegmentList( iter = ( (x,  x + int(numpy.random.randn() * 10.0 + 100.0)  ) \
                                                 for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )
        
        sampler = gat.SamplerAnnotator()
        
        counts = numpy.zeros( 1000 * nsegments, numpy.int )

        for x in range( self.ntests):
            sample = sampler.sample( segments, workspace )
            self.assertEqual( sample.sum(), segments.sum() )
            for start, end in sample: counts[start:end] += 1

        counts_within_workspace = []
        for start, end in workspace:
            counts_within_workspace.extend( counts[start:end] )

        # expected density: ntests * segment_size / workspace_size = numtest / 100
        self.assertAlmostEqual( numpy.mean(counts_within_workspace), 
                               self.ntests * segments.sum() / workspace.sum(),
                               places = 0 )
        print "standard deviation", numpy.std( counts_within_workspace )

    def testLengthSampling( self ):
        
        nsegments = 10000
        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        segments = gat.SegmentList( iter = [ (x, x + numpy.random.randn() * 10.0 + 100.0  ) \
                                                 for x in range(0, 1000 * nsegments, 1000) ],
                                    normalize = True )

        histogram = segments.getLengthDistribution( 1, 1000 * nsegments )
        hs = gat.HistogramSampler( histogram, 1 )

        samples = [hs.sample() for x in range(0, 1 * nsegments )]

        self.assertAlmostEqual( numpy.mean(samples),
                                100.0,
                                places = 0 )

        self.assertAlmostEqual( numpy.std(samples),
                                10.0,
                                places = 0 )

    def testPositionSampling( self ):
        '''test if we sample the exactly right amount of nucleoutides.'''
        workspace = gat.SegmentList( iter = ( (x, x + 100 ) for x in range( 0, 10000, 1000) ),
                                    normalize = True )
        segments = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 0, 10000, 1000) ),
                                    normalize = True )
        
        sampler = gat.SamplerAnnotator()
        
        counts = numpy.zeros( 10000, numpy.int )

        for x in range( self.ntests):
            sample = sampler.sample( segments, workspace )
            self.assertEqual( sample.sum(), segments.sum() )
            for start, end in sample: counts[start:end] += 1

        counts_within_workspace = []
        for start, end in workspace:
            counts_within_workspace.extend( counts[start:end] )

        # expected density: ntests * segment_size / workspace_size = numtest / 100
        self.assertAlmostEqual( numpy.mean(counts_within_workspace), 
                                self.ntests * segments.sum() / workspace.sum(),
                                places = 2 )
        print numpy.std( counts_within_workspace )

if __name__ == '__main__':
    unittest.main()
