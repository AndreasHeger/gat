'''test gat correctness
'''

import unittest
import random, tempfile, shutil, os, re, gzip, sys, subprocess
import gat
import numpy, math
import matplotlib.pyplot as plt

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

class TestSamplerLength( unittest.TestCase ):
    '''sample from length distribution.'''
    ntests = 1000
    nsegments = 10000
    sampler = None

    def setUp( self ):
        pass

    def testSamplingNormalDistribution( self ):

        if not self.sampler: return

        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        self.segments = gat.SegmentList( iter = [ (x, x + numpy.random.randn() * 10.0 + 100.0  ) \
                                                      for x in range(0, 1000 * self.nsegments, 1000) ],
                                         normalize = True )

        self.histogram = self.segments.getLengthDistribution( 1, 1000 * self.nsegments )

        hs = self.sampler( self.histogram, 1 )

        samples = [hs.sample() for x in range(0, 1 * self.nsegments )]

        delta = 1.0 

        self.assert_( abs( numpy.mean(samples) - 100.0) < delta )
        self.assert_( abs( numpy.std(samples) - 10.0) < delta )

    def testSamplingSNPs( self ):

        if not self.sampler: return

        # create normaly distributed lengths of mean 100.0 and sigma = 10.0
        self.segments = gat.SegmentList( iter = [ (x, x + 1  ) \
                                                      for x in range(0, 1000 * self.nsegments, 1000) ],
                                         normalize = True )

        self.histogram = self.segments.getLengthDistribution( 1, 1000 * self.nsegments )

        hs = self.sampler( self.histogram, 1 )

        samples = [hs.sample() for x in range(0, 1 * self.nsegments )]

        self.assertAlmostEqual( numpy.mean(samples),
                                1.0,
                                places = 0 )

        self.assertAlmostEqual( numpy.std(samples),
                                0.0,
                                places = 0 )
    
class TestSamplerLengthFast( TestSamplerLength ):

    sampler = gat.HistogramSampler

class TestSamplerLengthSlow( TestSamplerLength ):

    sampler = gat.HistogramSamplerSlow


class TestSamplerPosition( unittest.TestCase ):
    '''test segment position sampling (fixed length).'''

    ntests = 100000

    def getWorkspaceCounts( self, 
                            workspace, 
                            sampler, 
                            sample_length ):

        l = workspace.max()
        counts = numpy.zeros( l, numpy.int )

        for x in range( self.ntests):
            start, end, overlap = sampler.sample( sample_length )
            start = max(0, start)
            end = min( end, l )
            self.assert_( overlap > 0 )
            counts[start:end] += 1

        counts_within_workspace = []
        for start, end in workspace:
            counts_within_workspace.extend( counts[start:end] )

        if l > 10:
            dx = 10
        else:
            dx = 1

        counts_within_workspace = numpy.array( counts_within_workspace, dtype = numpy.int )
        newy = smooth( counts_within_workspace, window_len = dx)

        plt.figure()
        plt.plot( xrange(len(counts_within_workspace)), counts_within_workspace, '.', label = "coverage" )

        plt.plot( xrange(len(counts_within_workspace)), newy, '-', 
                  label="smooth - window = %i" % dx )

        plt.title( "%s : nworkspaces=%i, sample_length=%i" % ( str(self), len(workspace), sample_length ) )
        plt.xlabel( "position" )
        plt.ylabel( "counts" )
        plt.legend()
        plt.savefig( "test_%s.png" % re.sub( "[ ()]", "", str(self) ))

        # can be simplified
        expected = self.ntests * sample_length / float( workspace.sum()) * \
            float(workspace.sum()) / (workspace.sum() + sample_length  * len(workspace) )

        d = abs(counts_within_workspace.mean() - expected) / float(expected)
        self.assert_( d < 0.1, "expected counts (%f) != sampled counts (%f" % (expected,
                                                                               counts_within_workspace.mean()))
        
        d = numpy.std( counts_within_workspace )

        return counts_within_workspace

    def testPositionSamplingSingleWorkspace( self ):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace = gat.SegmentList( iter = [ (0,10000) ],
                                     normalize = True )

        sampler = gat.SegmentListSampler( workspace )
        
        counts_within_workspace = self.getWorkspaceCounts( workspace, sampler, 10 )

    def testPositionSamplingSplitWorkspace( self ):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace = gat.SegmentList( iter = [ (0,1000), (9000,10000) ],
                                     normalize = True )

        sampler = gat.SegmentListSampler( workspace )
        
        counts_within_workspace = self.getWorkspaceCounts( workspace, sampler, 10 )

    def testPositionSamplingSplitWorkspace2( self ):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        # note if the the space is too small, there are problems
        # in counting the expected overlap.
        workspace = gat.SegmentList( iter = [ (0,100), (200,300) ],
                                     normalize = True )

        sampler = gat.SegmentListSampler( workspace )

        counts_within_workspace = self.getWorkspaceCounts( workspace, 
                                                           sampler, 
                                                           2 )

    def testMultipleWorkspaces( self ):
        '''test if we sample the exactly right amount of nucleoutides
        and bases are overlapped uniformly.
        '''

        workspace = gat.SegmentList( iter = ( (x, x + 100 ) for x in range( 0, 10000, 1000) ),
                                    normalize = True )

        sampler = gat.SegmentListSampler( workspace )

        counts_within_workspace = self.getWorkspaceCounts( workspace, sampler, 10 )        

    def testSNPPositionSampling( self ):
        '''test if we sample the exactly right amount of nucleoutides.'''

        workspace = gat.SegmentList( iter = ( (x, x + 100 ) for x in range( 0, 10000, 1000) ),
                                    normalize = True )

        sampler = gat.SegmentListSampler( workspace )

        counts_within_workspace = self.getWorkspaceCounts( workspace, sampler, 1 )        

def getWorkspaceCounts( workspace, 
                        samples,
                        filename ):
    '''compute sample counts within workspace.

    returns array of size workspace with counts
    '''

    l = workspace.max()
    counts = numpy.zeros( l, numpy.int )
    ntests = len(samples)

    segment_sizes = []

    for sample_id, s in enumerate(samples):
        for start, end in s:
            start = max(0, start)
            end = min( end, l )
            counts[start:end] += 1
        ss = [x[1] - x[0] for x in s ] 

        segment_sizes.append( (numpy.std(ss),
                               min(ss), 
                               max(ss),
                               numpy.mean(ss),
                               numpy.median(ss) ) )
                               

    counts_within_workspace = []
    for start, end in workspace:
        counts_within_workspace.extend( counts[start:end] )

    if l > 10:
        dx = 10
    else:
        dx = 1

    counts_within_workspace = numpy.array( counts_within_workspace, dtype = numpy.int )
    newy = smooth( counts_within_workspace, window_len = dx)

    plt.figure()
    plt.subplot( "211" )
    plt.plot( xrange(len(counts_within_workspace)), counts_within_workspace, '.', label = "coverage" )

    plt.plot( xrange(len(counts_within_workspace)), newy, '-', 
              label="smooth - window = %i" % dx )

    plt.title( "%s : nworkspaces=%i, sample_length=%i" % ( filename, len(workspace), len(samples) ) )
    plt.xlabel( "position" )
    plt.ylabel( "counts" )
    plt.legend()

    plt.subplot( "212" )
    segment_sizes.sort()
    segment_sizes = zip(*segment_sizes)
    plt.plot( segment_sizes[0], label="stddev" )
    plt.plot( segment_sizes[1], label="min" )
    plt.plot( segment_sizes[2], label="max" )
    plt.plot( segment_sizes[3], label="mean" )
    plt.plot( segment_sizes[4], label="median" )
    plt.legend()
    plt.xlabel( "sample" )
    plt.ylabel( "segment size" )
    plt.savefig( filename )

    return counts_within_workspace

class TestSegmentSamplingGat( unittest.TestCase ):

    ntests = 1000

    def setUp(self):
        self.sampler = gat.SamplerAnnotator()

    def getSample( self, segments, workspace ):
        '''return a single sample'''
        return self.sampler.sample( segments, workspace )

    def getSamples( self, segments, workspace ):
        '''return a list of samples.'''
        samples = []

        for x in xrange( self.ntests):
            sample = self.getSample( segments, workspace )
            samples.append( sample )

        return samples

    def checkSample( self, samples, segments, workspace ):
        '''check if sample corresponds to expectation.

        segment_length: total length of segments (nsegments * segment_length)
        '''
        
        filename = "%s.png" % re.sub( "[ ()]", "", str(self) )
        filename = re.sub( "__main__", "", filename )

        self.assertEqual( self.ntests, len(samples) )

        counts_within_workspace = getWorkspaceCounts( workspace, 
                                                      samples,
                                                      filename )
        
        # check if each sample has the correct number of nucleotides  
        sums = [x.sum() for x in samples ]

        for x, s in enumerate(sums):
            self.assertEqual( s, segments.sum(), 
                              "incorrect number of nucleotides in sample %i, got %i, expected %i" %\
                                  (x, s, segments.sum() ) )


        # check if coverage is uniform over workspace
        # expected coverage
        # number_of_samples * segment_length / workspace_length
        # modifier: 
        # can be simplified
        nsegments = float(len(segments))
        segment_length = segments.sum() / nsegments
        expected = len(samples) * segments.sum() / float( workspace.sum())
        #   float(workspace.sum()) / (workspace.sum() + segments.sum()  * len(workspace) )
         

        d = abs(counts_within_workspace.mean() - expected) / float(expected)

        self.assert_( d < 0.1, "expected counts (%f) != sampled counts (%f)" % (expected,
                                                                                counts_within_workspace.mean()))
        
        stddev = numpy.std( counts_within_workspace )
        d = stddev / float(expected)

        self.assert_( d < 0.1, "coverage variation too large : stddev (%f) / %f = %f > 0.01" %\
                          ( d, 
                            stddev,
                            expected) )
        

    def testSegmentedWorkspaceSmallGap( self ):
        '''test sampling within a segmented workspace.

        The gap between adjacent workspaces is a single nucleotide
        to avoid the merging of intervals.
        '''
        nsegments = 10
        segment_size = 100

        workspace = gat.SegmentList( iter = ( (x, x + 999 ) \
                                                  for x in range( 0, 1000 * nsegments, 1000) ), 
                                     normalize = True )

        segments = gat.SegmentList( iter = ( (x,  x + segment_size)   \
                                                 for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )
        
        samples = self.getSamples( segments, workspace )

        self.checkSample( samples, segments, workspace )

    def testSegmentedWorkspaceLargeGap( self ):
        '''test sampling within a segmented workspace.

        The gap is the size of a segment.
        '''
        nsegments = 10
        segment_size = 100
        workspace = gat.SegmentList( iter = ( (x, x + 900 ) \
                                                  for x in range( 0, 1000 * nsegments, 1000) ), 
                                     normalize = True )

        segments = gat.SegmentList( iter = ( (x,  x + segment_size)   \
                                                 for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )
        
        samples = self.getSamples( segments, workspace )

        self.checkSample( samples, segments, workspace )

    def testSingleWorkspace( self ):
        '''test sampling within a single continuous workspace.'''
        nsegments = 10
        segment_size = 100

        workspace = gat.SegmentList( iter = [ (0, 10000) ],
                                     normalize = True )

        segments = gat.SegmentList( iter = ( (x,  x + segment_size)   \
                                                 for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )
        
        samples = self.getSamples( segments, workspace )

        self.checkSample( samples, segments, workspace )

    def testSegmentedWorkspace2x( self ):
        '''test if workspace segments are only twice the size of segmments.
        '''

        nsegments = 10 
        segment_size = 100

        workspace = gat.SegmentList( iter = ( (x, x + 2 * segment_size )
                                              for x in range( 0, 1000 * nsegments, 1000) ), 
                                     normalize = True )
        segments = gat.SegmentList( iter = ( (x, x + segment_size ) 
                                             for x in range( 0, 1000 * nsegments, 1000) ),
                                    normalize = True )

        samples = self.getSamples( segments, workspace )

        self.checkSample( samples, segments, workspace )

    # def testVariableLengthSampling( self ):
        
    #     nsegments = 10000
    #     # create normaly distributed lengths of mean 100.0 and sigma = 10.0
    #     segments = gat.SegmentList( iter = [ (x, x + numpy.random.randn() * 10.0 + 100.0  ) \
    #                                              for x in range(0, 1000 * nsegments, 1000) ],
    #                                 normalize = True )

    #     histogram = segments.getLengthDistribution( 1, 1000 * nsegments )
    #     hs = gat.HistogramSampler( histogram, 1 )

    #     samples = [hs.sample() for x in range(0, 1 * nsegments )]

    #     self.assertAlmostEqual( numpy.mean(samples),
    #                             100.0,
    #                             places = 0 )

    #     self.assertAlmostEqual( numpy.std(samples),
    #                             10.0,
    #                             places = 0 )


class TestSegmentSamplingTheAnnotator( TestSegmentSamplingGat ):
    '''use annotator to sample segments.'''
    
    cmd = '''/cpp-software/bin/java -Xmx8000M -cp /home/andreas/projects/annotator/lib/commons-cli-1.0.jar:/home/andreas/projects/annotator/lib/Annotator.jar app.Annotator -verbose 4'''

    def writeSegments( self, outfile, segmentlist, section ):
        
        if section == "workspace":
            prefix = "##Work"
        elif section == "segments":
            prefix = "##Seg"
        elif section == "annotations":
            prefix = "##Id\t0"

        contig = "chr1" 
        os.write( outfile,
                  "%s\t%s\t%s\n" % (
                prefix,
                contig,
                "\t".join( ["(%i,%i)" % (x,y) for x,y in segmentlist] ) ))

        if section == "annotations":
            os.write(outfile, "##Ann\ttest\t0" )

        os.close(outfile)

    def getSamples( self, 
                    segments,
                    workspace ):

        fsegments, nsegments = tempfile.mkstemp()
        fworkspace, nworkspace = tempfile.mkstemp()
        fannotations, nannotations = tempfile.mkstemp()

        self.writeSegments( fsegments, segments, "segments" )
        self.writeSegments( fworkspace, workspace, "workspace" )
        self.writeSegments( fannotations, workspace, "annotations" )

        statement = " ".join( (self.cmd,
                               "-dumpsegments",
                               "-iterations %i" % self.ntests,
                               "-annotation %s" % nannotations,
                               "-segments %s" % nsegments,
                               "-workspace %s" % nworkspace ) )

        process = subprocess.Popen(  statement,
                                     cwd = os.getcwd(), 
                                     shell = True,
                                     stdin = subprocess.PIPE,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE )

        # process.stdin.close()
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))

        samples = []
        for line in stdout.split("\n"):
            if not line.startswith( "##Segments" ): continue
            data = line.split("\t")[2:]
            coords = [ map(int, x[1:-1].split(",")) for x in data ]
            samples.append( gat.SegmentList( iter = coords ) )

        os.unlink( nsegments )
        os.unlink( nworkspace )
        os.unlink( nannotations )

        return samples

class TestStatsSNPSampling( unittest.TestCase ):
    '''test Stats by running a SNP (1-size interval) analysis.

    For SNPs, the hypergeometric distribution applies.
    '''

    sample_size = 10

    def check( self, workspace, annotations, segments ):

        workspace_size = workspace["chr1"].sum()
        
        sampler = gat.SamplerAnnotator( bucket_size = 1, nbuckets = workspace_size )

        counter = gat.CounterNucleotideOverlap()

        #print segments["default"]["chr1"]
        #print workspace["chr1"]
    
        annotator_results = gat.run( segments,
                                     annotations,
                                     workspace,
                                     sampler,
                                     counter,
                                     num_samples = self.sample_size )
        

        outfile = sys.stdout

        self.assertEqual( workspace_size, workspace["chr1"].sum() )

        tmp = gat.SegmentList( clone = segments["default"]["chr1"] )
        tmp.intersect( workspace["chr1"] )
        segment_size = tmp.sum()

        anno_mean, anno_pvalue, anno_std = [], [], []
        dist_mean, dist_pvalue, dist_std = [], [], []

        for track, r in annotator_results.iteritems():
            for annotation, result in r.iteritems():
                #print workspace["chr1"]
                #print annotations[annotation]["chr1"]
                
                intersection = gat.SegmentList( clone = annotations[annotation]["chr1"] )
                intersection.intersect( workspace["chr1"] )

                annotation_size = intersection.sum()

                # expected overlap for sampling with replacement
                expected_with = annotation_size / float(workspace_size)
                
                if annotation_size < workspace_size:
                    
                    # test sampling (# of expected)
                    # sampling without replacement follows hypergeometric distribution
                    # good = annotation_size
                    # bad = workspace - annotation_size
                    hyper = numpy.random.hypergeometric(annotation_size, 
                                                        workspace_size - annotation_size,
                                                        segment_size, 
                                                        self.sample_size )
                    hyper.sort()

                    # m = annotation_size # (white balls)
                    # N = workspace_size # (total balls)
                    # n = segment_size # (balls taken)
                    # variance = float(n * m * ( N - n ) * ( N -m )) / (N * N * (N - 1 )  )

                    # expected overlap for sampling without replacement
                    expected_without = hyper.mean()
                    error = hyper.std() * 4 # / segment_size

                    expected_std = hyper.std()

                    expected_pvalue = gat.getNPTwoSidedPValue( hyper, result.observed )

                    # print "\t".join( map(str, (result, 
                    #                            expected_without,
                    #                            expected_std,
                    #                            expected_pvalue,
                    #                            workspace_size,
                    #                            annotation_size) ) )

                    # for small sample size there might be no positive samples
                    if error == 0 and segment_size < 3: continue
                else:
                    # annotations = workspace
                    expected_without = segment_size 
                    expected_std = 0
                    expected_error = 0
                    expected_pvalue = 1.0
                    error = 0.1

                self.assert_( abs( result.expected - expected_without) < error,
                              "simulated results deviates from hypergeometric expectation for annotation `%s`: sizes(seg=%i/anno=%i/work=%i) observed=%f, expected=%f (%f, margin=%f)" %\
                                  ( annotation,
                                    segment_size,
                                    annotation_size,
                                    workspace_size,
                                    result.expected, expected_without, 
                                    result.expected - expected_without,
                                    error) )

                anno_mean.append( result.expected )
                anno_std.append( result.stddev )
                anno_pvalue.append( result.pvalue )

                dist_mean.append( expected_without )
                dist_std.append( expected_std )
                dist_pvalue.append( expected_pvalue )
                
                # plt.figure()
                # hist1, bins1 = numpy.histogram( hyper, new = True, bins=xrange(0, segment_size+10 ))
                # hist2, bins2 = numpy.histogram( result.samples, new = True, bins=bins1 )
                # plt.title( "%i - %i - %i" % (segment_size, 
                #                              annotation_size,
                #                              workspace_size) )
                # plt.plot( bins1[:-1], hist1, label = "dist" )
                # plt.plot( bins2[:-1], hist2, label = "simulated" )
                # plt.legend()
                # plt.show()

        plt.figure()
        plt.subplot( 221 )
        plt.scatter( anno_mean, dist_mean )
        plt.xlabel( "simulated - mean" )
        plt.ylabel( "distribution - mean" )
        plt.plot( anno_mean, anno_mean, "b" )
        plt.subplot( 222 )
        plt.scatter( anno_std, dist_std)
        plt.xlabel( "simulated - std" )
        plt.ylabel( "distribution - std" )
        plt.plot( anno_std, anno_std, "b" )
        plt.subplot( 223 )
        plt.scatter( anno_pvalue, dist_pvalue)
        plt.xlabel( "simulated - pvalue" )
        plt.ylabel( "distribution - pvalue" )
        plt.plot( anno_pvalue, anno_pvalue, "b" )
        plt.savefig( "test_%s.png" % re.sub( "[ ()]", "", str(self) ))

    def testSingleSNP( self ):

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
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
            
        self.check( workspace, annotations, segments )

    def testMultipleSNPsFullOverlap( self ):

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # 10 snps
        segments.add( "default", "chr1", gat.SegmentList( iter = [(x,x+1) for x in range(0,10)],
                                                          normalize = True ) )

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(10, 110, 5 ):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(0,y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testMultipleSNPsPartialOverlap( self ):
        '''test with multiple snps and decreasing
        amount of overlap with annotations.

        Tests if p-values are computed correctly.
        '''
        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        nsnps = 100

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # 10 snps
        segments.add( "default", "chr1", gat.SegmentList( iter = [(x,x+1) for x in range(0,nsnps)],
                                                          normalize = True ) )

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, nsnps, 1 ):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(y,nsnps+y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testIntervalsPartialOverlap( self ): 
        '''test with intervals with 
        increasing amount of overlap.

        '''
        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        size = 100

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # segment of size 10
        segments.add( "default", "chr1", gat.SegmentList( iter = [(0,size), ],
                                                          normalize = True ))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, size):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(y,size+y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testWorkspaces( self ):
        '''
        input:
            workspace = 500 segments of size 1000, separated by a gap of 1000
            annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
            segments = a SNP every 100 bp
           
        output:
        '''

        workspace_size = 1000000

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        # workspace of size 1000000
        workspaces.add( "default", "chr1", 
                        gat.SegmentList( iter = [(x, x+1000) for x in xrange(0,workspace_size, 2000)],
                                         normalize = True ) )

        # SNPs every 100bp
        segments.add( "default", "chr1", 
                      gat.SegmentList( iter = [(x,x+1) for x in xrange(0, workspace_size, 100) ],
                                       normalize = True ))
        
        annotations.add( "%03i" % start, "chr1",
                         gat.SegmentList( iter = [(0,workspace_size)],
                                          normalize = True ) ) 
            
        self.check( workspaces["default"], annotations, segments )

    def testFullAnnotation( self ):
        
        workspace_size = 1000000

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        # workspace of size 1000000
        workspaces.add( "default", "chr1", 
                        gat.SegmentList( iter = [(x, x+1000) for x in xrange(0,workspace_size, 2000)],
                                         normalize = True ) )

        # SNPs every 100bp
        segments.add( "default", "chr1", 
                      gat.SegmentList( iter = [(x,x+1) for x in xrange(0, workspace_size, 100) ],
                                       normalize = True ))
        
        size = 1000
        for start in xrange(0, size, 100):
            annotations.add( "%03i" % start, "chr1",
                             gat.SegmentList( iter = [(start+x,start+x+size) \
                                                          for x in xrange( 0, workspace_size, 2000)],
                                              normalize = True ) ) 
            
        self.check( workspaces["default"], annotations, segments )
        
class TestCaching( unittest.TestCase ):

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


class TestStats( unittest.TestCase ):

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

class TestStatsSampling( unittest.TestCase ):
    '''test Stats by running a SNP (1-size interval) analysis.

    For SNPs, the hypergeometric distribution applies.
    '''

    sample_size = 10

    def getSamples( self, segments,
             annotations,
             workspace,
             sampler,
             counter ):
        return gat.run( segments,
                        annotations,
                        workspace,
                        sampler,
                        counter,
                        num_samples = self.sample_size )
        
    def check( self, workspace, annotations, segments ):

        workspace_size = workspace["chr1"].sum()
        
        sampler = gat.SamplerAnnotator( bucket_size = 1, nbuckets = workspace_size )

        counter = gat.CounterNucleotideOverlap()

        #print segments["default"]["chr1"]
        #print workspace["chr1"]
    
        annotator_results = self.getSamples( segments,
                                             annotations,
                                             workspace,
                                             sampler,
                                             counter )

        outfile = sys.stdout

        self.assertEqual( workspace_size, workspace["chr1"].sum() )

        tmp = gat.SegmentList( clone = segments["default"]["chr1"] )
        tmp.intersect( workspace["chr1"] )
        segment_size = tmp.sum()

        anno_mean, anno_pvalue, anno_std = [], [], []
        dist_mean, dist_pvalue, dist_std = [], [], []

        for track, r in annotator_results.iteritems():
            for annotation, result in r.iteritems():
                #print workspace["chr1"]
                #print annotations[annotation]["chr1"]
                
                intersection = gat.SegmentList( clone = annotations[annotation]["chr1"] )
                intersection.intersect( workspace["chr1"] )

                annotation_size = intersection.sum()

                # expected overlap for sampling with replacement
                expected_with = annotation_size / float(workspace_size)
                
                if annotation_size < workspace_size:
                    
                    # test sampling (# of expected)
                    # sampling without replacement follows hypergeometric distribution
                    # good = annotation_size
                    # bad = workspace - annotation_size
                    hyper = numpy.random.hypergeometric(annotation_size, 
                                                        workspace_size - annotation_size,
                                                        segment_size, 
                                                        self.sample_size )
                    hyper.sort()

                    # m = annotation_size # (white balls)
                    # N = workspace_size # (total balls)
                    # n = segment_size # (balls taken)
                    # variance = float(n * m * ( N - n ) * ( N -m )) / (N * N * (N - 1 )  )

                    # expected overlap for sampling without replacement
                    expected_without = hyper.mean()
                    error = hyper.std() * 4 # / segment_size

                    expected_std = hyper.std()

                    expected_pvalue = gat.getNPTwoSidedPValue( hyper, result.observed )

                    # print "\t".join( map(str, (result, 
                    #                            expected_without,
                    #                            expected_std,
                    #                            expected_pvalue,
                    #                            workspace_size,
                    #                            annotation_size) ) )

                    # for small sample size there might be no positive samples
                    if error == 0 and segment_size < 3: continue
                else:
                    # annotations = workspace
                    expected_without = segment_size 
                    expected_std = 0
                    expected_error = 0
                    expected_pvalue = 1.0
                    error = 0.1

                self.assert_( abs( result.expected - expected_without) < error,
                              "simulated results deviates from hypergeometric expectation for annotation `%s`: sizes(seg=%i/anno=%i/work=%i) observed=%f, expected=%f (%f, margin=%f)" %\
                                  ( annotation,
                                    segment_size,
                                    annotation_size,
                                    workspace_size,
                                    result.expected, expected_without, 
                                    result.expected - expected_without,
                                    error) )

                anno_mean.append( result.expected )
                anno_std.append( result.stddev )
                anno_pvalue.append( result.pvalue )

                dist_mean.append( expected_without )
                dist_std.append( expected_std )
                dist_pvalue.append( expected_pvalue )
                
                # plt.figure()
                # hist1, bins1 = numpy.histogram( hyper, new = True, bins=xrange(0, segment_size+10 ))
                # hist2, bins2 = numpy.histogram( result.samples, new = True, bins=bins1 )
                # plt.title( "%i - %i - %i" % (segment_size, 
                #                              annotation_size,
                #                              workspace_size) )
                # plt.plot( bins1[:-1], hist1, label = "dist" )
                # plt.plot( bins2[:-1], hist2, label = "simulated" )
                # plt.legend()
                # plt.show()

        plt.figure()
        plt.subplot( 221 )
        plt.scatter( anno_mean, dist_mean )
        plt.xlabel( "simulated - mean" )
        plt.ylabel( "distribution - mean" )
        plt.plot( anno_mean, anno_mean, "b" )
        plt.subplot( 222 )
        plt.scatter( anno_std, dist_std)
        plt.xlabel( "simulated - std" )
        plt.ylabel( "distribution - std" )
        plt.plot( anno_std, anno_std, "b" )
        plt.subplot( 223 )
        plt.scatter( anno_pvalue, dist_pvalue)
        plt.xlabel( "simulated - pvalue" )
        plt.ylabel( "distribution - pvalue" )
        plt.plot( anno_pvalue, anno_pvalue, "b" )
        plt.savefig( "test_%s.png" % re.sub( "[ ()]", "", str(self) ))

    def testSingleSNP( self ):

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
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
            
        self.check( workspace, annotations, segments )

    def testMultipleSNPsFullOverlap( self ):

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # 10 snps
        segments.add( "default", "chr1", gat.SegmentList( iter = [(x,x+1) for x in range(0,10)],
                                                          normalize = True ) )

        # annotations: a collection of segments with increasing density
        # all are overlapping the segments
        for y in range(10, 110, 5 ):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(0,y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testMultipleSNPsPartialOverlap( self ):
        '''test with multiple snps and decreasing
        amount of overlap with annotations.

        Tests if p-values are computed correctly.
        '''
        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        nsnps = 100

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # 10 snps
        segments.add( "default", "chr1", gat.SegmentList( iter = [(x,x+1) for x in range(0,nsnps)],
                                                          normalize = True ) )

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, nsnps, 1 ):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(y,nsnps+y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testIntervalsPartialOverlap( self ): 
        '''test with intervals with 
        increasing amount of overlap.

        '''
        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        workspace_size = 1000

        size = 100

        # workspace of size 1000000
        workspaces.add( "default", "chr1", gat.SegmentList( iter = [(0,workspace_size),],
                                                            normalize = True ) )
        workspace = workspaces["default"]

        # segment of size 10
        segments.add( "default", "chr1", gat.SegmentList( iter = [(0,size), ],
                                                          normalize = True ))

        # annotations: a collection of segments.
        # overlap increases
        for y in range(0, size):
            annotations.add( "%03i" % y, "chr1",
                             gat.SegmentList( iter = [(y,size+y),],
                                              normalize = True ) ) 
            
        self.check( workspace, annotations, segments )

    def testWorkspaces( self ):
        '''
        input:
            workspace = 500 segments of size 1000, separated by a gap of 1000
            annotations = 500 segments of size 1000, separated by a gap of 1000, shifted up 100 bases
            segments = a SNP every 100 bp
           
        output:
        '''

        workspace_size = 1000000

        workspaces, segments, annotations = \
            gat.IntervalCollection( "workspace" ), \
            gat.IntervalCollection( "segment" ), \
            gat.IntervalCollection( "annotation" )

        # workspace of size 1000000
        workspaces.add( "default", "chr1", 
                        gat.SegmentList( iter = [(x, x+1000) for x in xrange(0,workspace_size, 2000)],
                                         normalize = True ) )

        # SNPs every 100bp
        segments.add( "default", "chr1", 
                      gat.SegmentList( iter = [(x,x+1) for x in xrange(0, workspace_size, 100) ],
                                       normalize = True ))
        
        size = 1000
        for start in xrange(0, size, 100):
            annotations.add( "%03i" % start, "chr1",
                             gat.SegmentList( iter = [(start+x,start+x+size) \
                                                          for x in xrange( 0, workspace_size, 2000)],
                                              normalize = True ) ) 
            
        self.check( workspaces["default"], annotations, segments )


if __name__ == '__main__':
    unittest.main()
