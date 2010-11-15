"""test gat implementation."""

import unittest
import random, tempfile, shutil, os, re, gzip, sys
import gat
import numpy, math

class TestSegmentList( unittest.TestCase ):
    '''test segment list implementation.'''

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

    def testNormalize1b( self ):
        '''non-overlapping segments.'''

        ss = [ (x, x + 10 ) for x in range( 100, 1100, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()
        print str(s)
        self.assertEqual( len(s), 10 )
        self.assertEqual( s.sum(), 100 )

        

    def testNormalizeEmpty( self ):
        '''non-overlapping segments.'''

        s = gat.SegmentList()
        self.assertEqual( len(s), 0)
        s.normalize()
        self.assertEqual( len(s), 0)
        self.assertEqual( s.isNormalized, 1)

    def testNormalizeEmptySegment( self ):
        s = gat.SegmentList( iter = [(0, 0),] )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 0)

        s = gat.SegmentList( iter = [(0, 0),(0,0)] )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 0)

        s = gat.SegmentList( iter = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9)] )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 1)

    def testNormalize2( self ):
        '''overlapping segments.'''

        ss = [ (x, x + 1000 ) for x in range( 0, 1000, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()
        self.assertEqual( len(s), 1 )
        self.assertEqual( s.sum(), 1900 )

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

    def testNormalize5( self ):
        ss = [(489, 589), (1966, 2066), (2786, 2886), (0, 0), (3889, 3972), (3998, 4098), (6441, 6541), (6937, 7054), (7392, 7492), (8154, 8254), (9046, 9146)]
        s = gat.SegmentList( iter = ss )
        s.normalize()
        self.assertEqual( len(s), len( [x for x in ss if x[1]-x[0] > 0] ) )
        self.assertEqual( s.sum(), 1000 )

    def testExtend( self ):
        
        s1 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 0, 1000, 100) ] )
        s2 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 2000, 3000, 100) ] )
        s1.extend(s2 )
        self.assertEqual( s1.sum(), s2.sum() * 2 )
        self.assertEqual( len(s1), len(s2) * 2 )

    def testTrim( self ):
        '''test trimming over full range of insertion points and deletions.'''

        for point in xrange( 0, 1100 ):
            for size in xrange( 0, 200 ):
                s = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 0, 1000, 100) ],
                                     normalize = True )
                orig = s.sum()
                s.trim( point, size )
                self.assertEqual( orig - size, s.sum() ) 

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

    def testOverlapOutOfRange( self ):
        self.assertEqual( self.a.overlapWithRange( -100,5 ), 5)
        self.assertEqual( self.a.overlapWithRange( -100,-50 ), 0)
        self.assertEqual( self.a.overlapWithRange( 905,1100 ), 5)
        self.assertEqual( self.a.overlapWithRange( 1000,1100 ), 0)

    def testOverlapAll( self ):
        for x in range( 0, 1000, 100):
            for y in range( 0, 10 ):
                self.assertEqual( self.a.overlapWithRange( x+y,x+y+1), 1, \
                                      "no overlap failure at %i: %i" % (x+y, self.a.overlapWithRange( x+y,x+y+1)))
            for y in range( 10, 100 ):
                self.assertEqual( self.a.overlapWithRange( x+y,x+y+1), 0, "overlap failure at %i" % (x+y) )


class TestSegmentListIntersection( unittest.TestCase):

    def setUp( self ):
        #[(0, 10), (100, 110), (200, 210), (300, 310), (400, 410), (500, 510), (600, 610), (700, 710), (800, 810), (900, 910)]
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
        self.assertEqual( r.isEmpty, True )

    def testPartialIntersection( self ):
        b = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 5, 1000, 100) ), normalize = True )
        r = b.intersect( self.a )
        self.assertEqual( len(r), len(self.a) )
        self.assertEqual( r.sum(), self.a.sum() / 2 )

    def testOverlap( self ):
        '''test if number of segments intersection is correct.'''

        b = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 5, 1000, 100) ), normalize = True )
        self.assertEqual( self.a.intersectionWithSegments( b ), len(b) )
        self.assertEqual( b.intersectionWithSegments( self.a ), len(b) )

        # no intersection
        b = gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 10, 1000, 100) ), normalize = True )
        self.assertEqual( self.a.intersectionWithSegments( b ), 0 )
        self.assertEqual( b.intersectionWithSegments( self.a ), 0 )
        
        # double the number of segments in b
        b = gat.SegmentList( iter = [(x, x + 5 ) for x in range( 0, 1000, 100) ] +\
                                 [(x+5, x + 10 ) for x in range( 0, 1000, 100) ], \
                                 normalize = True )
        self.assertEqual( self.a.intersectionWithSegments( b ), 10 )
        self.assertEqual( b.intersectionWithSegments( self.a ), 20 )

    def testFilter( self ):
        
        b = gat.SegmentList( iter = ( (x, x + 5 ) for x in range( 500, 2000, 100) ), normalize = True )
        self.assertEqual( b.filter(self.a).asList(), [(500, 505), (600, 605), (700, 705), (800, 805), (900, 905)] )

class TestIntervalCollection( unittest.TestCase):

    def setUp( self ):
        self.a = gat.IntervalCollection( "a" )
        self.a.add( "track1", "contig1",
                    gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 0, 1000, 100) ), normalize = True ) )
        self.a.add( "track2", "contig1",
                    gat.SegmentList( iter = ( (x, x + 10 ) for x in range( 1000, 2000, 100) ), normalize = True ) )
        
    def testSaveLoad( self ):
        
        fn = "tmp_testSaveLoad.bed"
        outfile = open(fn, "w")
        self.a.save( outfile )
        outfile.close()

        b = gat.IntervalCollection("b")
        b.load( fn )

        self.assertEqual( len(self.a), len(b) )
        self.assertEqual( self.a.tracks, b.tracks )
        self.assertEqual( self.a.sum(), b.sum() )
        for t in b.tracks:
            self.assertEqual( sorted(self.a[t].keys()), sorted(b[t].keys()) )
            for x in self.a[t].keys():
                self.assertEqual( self.a[t][x], b[t][x] )

class TestSamples( unittest.TestCase ):
    
    nsamples = 1000
    ntracks = 50
    nsegments = 10000
    nisochores = 20


    def testDelete( self ):
        # from guppy import hpy
        # hp = hpy()
        # hp.setrelheap()

        samples = gat.Samples()

        for track in xrange(self.ntracks):
            track_id = str(track)
            print track_id
            for sample in xrange(self.nsamples):
                sample_id = str(sample)
                for isochore in xrange( self.nisochores):
                    isochore_id = str(isochore)
                    r = gat.SegmentList( allocate = self.nsegments )
                    samples.add( track_id, sample_id, isochore_id, r )

            del samples[track_id]
            print len(samples)
            # h = hp.heap()

            # print h


if __name__ == '__main__':
    unittest.main()
