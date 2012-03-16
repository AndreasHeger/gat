"""test gat implementation."""

import unittest
import random, tempfile, shutil, os, re, gzip, sys
import gat
import numpy, math

class GatTest( unittest.TestCase ):
    def shortDescription( self ): return None

class TestSegmentList( GatTest ):
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
        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testNormalize1b( self ):
        '''non-overlapping segments.'''

        ss = [ (x, x + 10 ) for x in range( 100, 1100, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()
        self.assertEqual( len(s), 10 )
        self.assertEqual( s.sum(), 100 )
        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testNormalizeEmpty( self ):
        '''non-overlapping segments.'''

        s = gat.SegmentList()
        self.assertEqual( len(s), 0)
        s.normalize()
        self.assertEqual( len(s), 0)
        self.assertEqual( s.isNormalized, 1)
        s2 = gat.SegmentList()
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testNormalizeEmptySegment( self ):
        s = gat.SegmentList( iter = [(0, 0),] )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 0)

        s = gat.SegmentList( iter = [(0, 0),(0,0)] )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 0)

        ss = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9)]
        s = gat.SegmentList( iter = ss )
        s.normalize()        
        self.assertEqual( s.isNormalized, 1)
        self.assertEqual( len(s), 1)

        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

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
        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testNormalize4( self ):
        # test multiple interleaved segments
        ss = [ (x, x + 100 ) for x in range( 0, 1000, 10) ]
        s = gat.SegmentList()
        for start, end in ss: s.add( start, end )
        s.normalize()
        self.assertEqual( len(s), 1 )
        self.assertEqual( s.sum(), 1090 )
        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testNormalize5( self ):
        ss = [(489, 589), (1966, 2066), (2786, 2886), (0, 0), (3889, 3972), (3998, 4098), (6441, 6541), (6937, 7054), (7392, 7492), (8154, 8254), (9046, 9146)]
        s = gat.SegmentList( iter = ss )
        s.normalize()
        self.assertEqual( len(s), len( [x for x in ss if x[1]-x[0] > 0] ) )
        self.assertEqual( s.sum(), 1000 )
        s2 = gat.SegmentList( iter = ss )
        s2.merge( distance = -1 )
        self.assertEqual( s, s2 )

    def testExtend( self ):
        
        s1 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 0, 1000, 100) ] )
        s2 = gat.SegmentList( iter =  [ (x, x + 100 ) for x in range( 2000, 3000, 100) ] )
        s1.extend(s2 )
        self.assertEqual( s1.sum(), s2.sum() * 2 )
        self.assertEqual( len(s1), len(s2) * 2 )

    def testTrim( self ):
        '''test trimming over full range of insertion points and deletions.'''

        for point in xrange( 0, 1000 ):
            for size in xrange( 0, 300 ):
                ss = [ (x, x + 100 ) for x in range( 0, 1000, 100) ]
                s = gat.SegmentList( iter = ss,
                                     normalize = True )
                orig = s.sum()
                s.trim( point, size )
                self.assertEqual( orig - size, s.sum(),
                                  "trimming error at %i:%i: expected %i, got %i, %s" %\
                                      (point, size,
                                       orig-size, s.sum(), str(s) ) )


    def testInsertionPoint( self ):
        '''check insertion point for normalized segment lists'''
        ss = [ (x, x + 10 ) for x in range( 0, 100, 10) ]
        s = gat.SegmentList( iter = ss,
                             normalize = True )
        
        for point in xrange( 0, 100 ):
            p = s.getInsertionPoint( point, point + 1 )
            self.assertEqual( p, point // 10 )

        ss = [ (x, x + 10 ) for x in range( 0, 100, 20) ]
        s = gat.SegmentList( iter = ss,
                             normalize = True )
        
        for point in xrange( 0, 100 ):
            p = s.getInsertionPoint( point, point + 1 )
            if point >= 90: 
                self.assertEqual( p, len(s) )
            else:
                self.assertEqual( p, point // 20 ) 

        ss = [ (x, x + 10 ) for x in range( 10, 100, 20) ]
        s = gat.SegmentList( iter = ss,
                             normalize = True )
        
        for point in xrange( 0, 100 ):
            p = s.getInsertionPoint( point, point + 1 )
            self.assertEqual( p, (point - 10) // 20 )

    def testInsertionPointNonNormalized( self ):
        '''check insertion point for unnormalized segment lists.'''
        ss = [ (x, x + 20 ) for x in range( 0, 100, 10) ]
        s = gat.SegmentList( iter = ss,
                             normalize = False )
        
        for point in xrange( 0, 100 ):
            p = s.getInsertionPoint( point, point + 1 )
            print point, p


    def testMergeAdjacent( self ):
        ss = [ (x, x + 100  ) for x in range( 0, 1000, 100) ]
        random.shuffle(ss)
        s = gat.SegmentList( iter = ss )
        s.merge( )
        self.assertEqual( len(s), 1 )
        self.assertEqual( s.sum(), 1000 )

    def testMergeNeighbours( self ):
  
        for y in range( 0,5 ):
            ss = [ (x, x + 100 - y ) for x in range( 0, 1000, 100) ]
            random.shuffle(ss)
            for x in range(0,y+1):
                s = gat.SegmentList( iter = ss )
                s.merge( distance = x )
                if x < y:
                    self.assertEqual( len(s), 10 )
                    self.assertEqual( s.sum(), 1000 - 10 * y)
                else:
                    self.assertEqual( len(s), 1 )
                    self.assertEqual( s.sum(), 1000 - y )

    def testExpand( self ):
        ss = [ (x, x + 10 ) for x in range( 0, 1000, 100) ]
        s = gat.SegmentList( iter = ss, normalize = True )
        s.expand_segments( 2.0 )
        self.assertEqual( s.sum(), 195)
        self.assertEqual( len(s), 10 )

    def testGetFilledSegmentsFromStart( self ):
        ss = [ (x, x + 10 ) for x in range( 0, 120, 20) ]
        s = gat.SegmentList( iter = ss, normalize = True )
        for x in range( 0, 120, 5 ):
            f = s.getFilledSegmentsFromStart( x, 20 )
            self.assertEqual( f.sum(), 20 )
            if x >= 110:
                self.assertEqual( f.min(), s.min() )
                self.assertEqual( f.max(), 30 )
            elif x > 80: 
                self.assertEqual( f.min(), s.min() )
                self.assertEqual( f.max(), s.max() )
            else:
                if x in (0,20,40,60,80):
                    self.assertEqual( f.min(), x )
                    self.assertEqual( f.max(), x + 30 )
                elif x in (10,30,50,70):
                    self.assertEqual( f.min(), x + 10)
                    self.assertEqual( f.max(), x + 40 )
                elif x in (5,25,45,65):
                    self.assertEqual( f.min(), x )
                    self.assertEqual( f.max(), x + 40 )
                elif x in (15,35,55,75):
                    self.assertEqual( f.min(), x + 5 )
                    self.assertEqual( f.max(), x + 35 )
            self.assertEqual( s.sum(), s.getFilledSegmentsFromStart(x, 100 ).sum() )

    def testGetFilledSegmentsFromEnd( self ):
        ss = [ (x, x + 10 ) for x in range( 0, 120, 20) ]
        s = gat.SegmentList( iter = ss, normalize = True )
        for x in range( 0, 120, 5 ):
            f = s.getFilledSegmentsFromEnd( x, 20 )
            self.assertEqual( f.sum(), 20 )
            if x == 0:
                self.assertEqual( f.max(), s.max() )
                self.assertEqual( f.min(), 80 )
            elif x < 30: 
                self.assertEqual( f.max(), s.max() )
                self.assertEqual( f.min(), s.min() )
            else:
                if x in (40,60,80,100):
                    self.assertEqual( f.min(), x - 40 )
                    self.assertEqual( f.max(), x - 10 )
                elif x in (30,50,70,90):
                    self.assertEqual( f.min(), x - 30 )
                    self.assertEqual( f.max(), x )
                elif x in (45,65,85):
                    self.assertEqual( f.min(), x - 40)
                    self.assertEqual( f.max(), x  )
                elif x in (35,55,75,95):
                    self.assertEqual( f.min(), x - 35 )
                    self.assertEqual( f.max(), x - 5 )
                    
            self.assertEqual( s.sum(), s.getFilledSegmentsFromEnd(x, 100 ).sum() )

class TestSegmentListOverlap( GatTest ):
    
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
        self.assertEqual( self.a.overlapWithRange( -100, 5 ), 5)
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


class TestSegmentListIntersection( GatTest):

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

        # test negative segments
        b = gat.SegmentList( iter = ( (-1,56), ) )
        c = gat.SegmentList( iter = [(0, 50), (75, 125)] )
        self.assertEqual( b.filter(c).asList(), [(-1,56)] )

        b = gat.SegmentList( iter = ( (-1,56), ) )
        c = gat.SegmentList( iter = [(-10, 10)] )
        self.assertEqual( b.filter(c).asList(), [(-1,56)] )

class TestIntervalCollection( GatTest):

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

        
class TestSamples( GatTest ):
    
    nsamples = 1000
    ntracks = 50
    nsegments = 10000
    nisochores = 20

    def testDelete( self ):
        '''test to track down a memory leak.'''
        # from guppy import hpy
        # hp = hpy()
        # hp.setrelheap()
        return

        samples = gat.Samples()

        for track in xrange(self.ntracks):
            track_id = str(track)
            # print track_id
            for sample in xrange(self.nsamples):
                sample_id = str(sample)
                for isochore in xrange( self.nisochores):
                    isochore_id = str(isochore)
                    r = gat.SegmentList( allocate = self.nsegments )
                    samples.add( track_id, sample_id, isochore_id, r )

            del samples[track_id]
            # print len(samples)
            # h = hp.heap()

            # print h


class TestCaching( GatTest ):

    sample_size = 10
    workspace_size = 1000

    def testCaching( self ):

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,self.workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        segments.add( "default", "chr1", gat.SegmentList( iter = [(0,1),],
                                                          normalize = True ) )

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(1, 100, 2 ):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(0,y),],
                                              normalize = True ) ) 
            
        workspace_size = workspace["chr1"].sum()

        sampler = gat.SamplerAnnotator( bucket_size = 1, nbuckets = self.workspace_size )
        
        if os.path.exists( "test.cache"): os.remove( "test.cache" )
        
        outsamples = gat.SamplesCached( "test.cache" )
        saved_samples = {}

        for track in segments.tracks:
            segs = segments[track]
            for x in xrange( self.sample_size ):
                for isochore in segs.keys():
                    r = sampler.sample( segs[isochore], workspace[isochore] )
                    saved_samples[(track,x,isochore)] = r
                    outsamples.add( track, x, isochore, r )
        
        del outsamples
        
        insamples = gat.SamplesCached( "test.cache" )
        
        for track in segments.tracks:
            segs = segments[track]
            for x in xrange( self.sample_size ):
                for isochore in segs.keys():
                    insamples.load( track, x, isochore )
                    
        # compare
        for track in segments.tracks:
            segs = segments[track]
            for x in xrange( self.sample_size ):
                for isochore in segs.keys():
                    self.assertEqual( saved_samples[(track,x,isochore)].asList(),
                                      insamples[track][x][isochore].asList() )


class TestStats( GatTest ):

    ntracks = 10 # 17
    nannotations = 10 # 90
    nsamples = 1000

    def testPValueComputation( self ):
        
        annotation_size = 100
        workspace_size = 1000
        segment_size = 10

        l = 10

        for y in xrange(1, l):

            samples = [1] * y + [0] * (l - y)

            for x, s in enumerate(samples):

                g = gat.AnnotatorResult( "track", "samples",
                                         s,
                                         samples )
                self.assertEqual( g.isSampleSignificantAtPvalue( x, g.pvalue ), True )
                                  

                t = 0
                for z, s2 in enumerate(samples):
                    t += g.isSampleSignificantAtPvalue( z, g.pvalue )
                fpr = float(t) / l

                # == should work, but there is a problem 
                # for pvalue = 0.5
                self.assert_( fpr >= g.pvalue - 0.0001, 
                              "fpr (%f) != pvalue (%f): y=%i, x=%i, s=%i, t=%i" % \
                                  ( fpr, g.pvalue, y, x, s, t) )


    def testPValueComputation2( self ):
        
        samples = [0] * 66 + [1] * 2 + [2] * 20 + [3] * 1 + [4] * 6 + [6] * 2 + [8] * 2 + [16] * 1
        
        obs = 16

        g = gat.AnnotatorResult( "track", "samples",
                                 obs,
                                 samples )

        self.assertEqual( g.pvalue, 0.01 )


    def testStats( self ):
        
        annotation_size = 100
        workspace_size = 1000
        segment_size = 10

        observed = numpy.random.hypergeometric(annotation_size, 
                                               workspace_size - annotation_size,
                                               segment_size, 
                                               self.ntracks * self.nannotations )

        results = []
        x = 0
        for track in range(self.ntracks ):
            for annotation in range(self.nannotations ):
                samples = numpy.random.hypergeometric(annotation_size, 
                                                      workspace_size - annotation_size,
                                                      segment_size, 
                                                      self.nsamples )

                samples.sort()
                # print "obs", observed[x], "sample=", samples

                results.append( gat.AnnotatorResult( str(track),                                                      str(annotation),
                                                     observed[x],
                                                     samples ) )

                x += 1

        gat.computeFDR( results )
        for r in results: 
            self.assert_( r.qvalue > 0.5, "%f" % r.qvalue  )



if __name__ == '__main__':
    unittest.main()
