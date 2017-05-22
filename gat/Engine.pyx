# cython: embedsignature=True
# cython: profile=False

import collections
import re
import os
import math
import random

import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats

cimport cython
from libc.stdio cimport FILE, fopen, fclose, feof
from libc.stdio cimport fread, fwrite, ftell, fseek, SEEK_SET
from libc.stdlib cimport realloc, malloc, calloc, free, atol
from libc.string cimport memcpy, memmove, memchr, strlen
from libc.math cimport floor
from libc.errno cimport errno
from posix.types cimport off_t
from posix.mman cimport mmap, munmap, shm_open, shm_unlink
from posix.mman cimport MAP_SHARED, PROT_READ, PROT_WRITE
from posix.stat cimport S_IRUSR, S_IWUSR
from posix.fcntl cimport O_CREAT, O_RDWR, O_RDONLY
from posix.unistd cimport ftruncate

# import SegmentList and PositionList
from CoordinateList cimport CoordinateList

cimport SegmentList
import SegmentList
from SegmentList cimport SegmentList, PositionDifference, Segment, Position, force_bytes, force_str

from PositionList cimport PositionList

# Numpy import
# note both import and cimport are necessary
import numpy
cimport numpy
DTYPE_INT = numpy.int
ctypedef numpy.int_t DTYPE_INT_t
DTYPE_FLOAT = numpy.float
ctypedef numpy.float_t DTYPE_FLOAT_t

#####################################################
## scipy import
## required for pvalues based on parametric distributions
try:
    import scipy
    import scipy.stats
    HAS_SCIPY = True
except:
    HAS_SCIPY = False

# min/max are not optimized, so declare them as C functions
# declare as signed comparisons as a Position might be negative
@cython.profile(False)
cdef inline PositionDifference lmin(PositionDifference a, PositionDifference b) nogil:
    if a < b:
        return a
    return b

@cython.profile(False)
cdef inline PositionDifference lmax(PositionDifference a, PositionDifference b) nogil:
    if a > b:
        return a
    return b

@cython.profile(False)
cdef inline double dmin(double a, double b) nogil:
    if a < b:
        return a
    return b

@cython.profile(False)
cdef inline double dmax(double a, double b) nogil:
    if a > b:
        return a
    return b

@cython.profile(False)
cdef inline PositionDifference segment_overlap( Segment a, Segment b) nogil:
    return lmax(0, <PositionDifference>lmin( a.end, b.end) - \
                <PositionDifference>lmax(a.start, b.start))

@cython.profile(False)
cdef inline PositionDifference range_overlap(
    Position astart, Position aend, Position bstart, Position bend) nogil:
    return lmax(0, 
                <PositionDifference>lmin( aend, bend) - 
                <PositionDifference>lmax(astart, bstart))

@cython.profile(False)
cdef inline PositionDifference segment_overlap_raw(Segment a, Segment b) nogil:
    return <PositionDifference>lmin(a.end, b.end) - \
        <PositionDifference>lmax(a.start, b.start)

@cython.profile(False)
cdef inline PositionDifference segment_length( Segment a) nogil:
    return <PositionDifference>a.end - <PositionDifference>a.start

# trick to permit const void * in function definitions
cdef extern from *:
    ctypedef void * const_void_ptr "const void*"

@cython.profile(False)
cdef int cmpSegments( const_void_ptr s1, const_void_ptr s2 ):
    return <PositionDifference>(<Segment *>s1).start - <PositionDifference>(<Segment *>s2).start

@cython.profile(False)
cdef int cmpSegmentsStartAndEnd( const_void_ptr s1, const_void_ptr s2 ):
    cdef int x
    x = <PositionDifference>(<Segment *>s1).start - <PositionDifference>(<Segment *>s2).start
    if x != 0: return x
    return <PositionDifference>(<Segment *>s1).end - <PositionDifference>(<Segment *>s2).end

@cython.profile(False)
cdef int cmpPosition( const_void_ptr s1, const_void_ptr s2 ):
    return (<Position*>s1)[0] - (<Position*>s2)[0]

@cython.profile(False)
cdef int cmpDouble( const_void_ptr s1, const_void_ptr s2 ):
    # see http://www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html
    cdef double da = (<double*>s1)[0]
    cdef double db = (<double*>s2)[0]
    return (da > db) - (da < db)

cdef class SegmentListSamplerSlow:

    cdef SegmentList segment_list
    cdef numpy.ndarray  cdf
    cdef Position total_size

    def __init__(self, SegmentList segment_list ):

        assert len(segment_list) > 0, "sampling from empty segment list"

        self.segment_list = segment_list
        self.cdf = numpy.cumsum( [x[1] - x[0] for x in self.segment_list ] )
        self.total_size = self.cdf[-1]

    cpdef sample(self, Position length):
        '''return a new position within segment list.'''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef Position r, offset, pos
        cdef PositionDifference overlap
        cdef size_t segment_index

        r = numpy.random.randint(0, self.total_size)
        segment_index = numpy.searchsorted(self.cdf, r)
        offset = r - self.cdf[segment_index]
        if offset == 0:
            pos = self.segment_list.segments[segment_index].start
        else:
            pos = self.segment_list.segments[segment_index].end + offset

        overlap = self.segment_list.overlapWithRange( pos, pos+length )

        assert overlap > 0, "sample %i does not overlap with workspace: offset=%i, r=%i, index=%i/%i" %\
            (pos, r - self.cdf[segment_index], r, segment_index, self.segment_list.nsegments)

        return pos, overlap

cdef class SegmentListSamplerWithEdgeEffects:

    cdef SegmentList segment_list
    cdef Position * cdf
    cdef Position total_size
    cdef Position nsegments

    def __init__(self, SegmentList segment_list):
        cdef Position i, totsize
        assert len(segment_list) > 0, "sampling from empty segment list"

        self.segment_list = segment_list
        self.nsegments = len(segment_list)
        self.cdf = <Position*>malloc(sizeof(Position) * self.nsegments)
        if not self.cdf:
            raise MemoryError(
                "out of memory when allocation %i bytes" %
                sizeof( sizeof(Position) * self.nsegments))
        self.total_size = 0
        for i from 0 <= i < len(segment_list):
            self.total_size += segment_length(segment_list.segments[i])
            self.cdf[i] = self.total_size

    cpdef sample(self, Position length):
        '''return a new position within segment list.

        This method both samples a position within the workspace
        and a position within length.
        '''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef Position rpos, rseg, offset, pos, start, end
        cdef PositionDifference overlap
        cdef size_t segment_index

        # r = rand() / (RAND_MAX / total_size + 1)
        rpos = numpy.random.randint( 0, self.total_size )
        segment_index = searchsorted( self.cdf,
                                      self.nsegments,
                                      sizeof(Position),
                                      &rpos,
                                      &cmpPosition,
                                      )

        offset = rpos - self.cdf[segment_index]
        if offset == 0:
            pos = self.segment_list.segments[segment_index+1].start
        else:
            pos = self.segment_list.segments[segment_index].end + offset

        rseg = numpy.random.randint(0, length)

        start, end = pos - rseg, pos - rseg + length

        # print "sample %i-%i, offset=%i, pos=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
        #     (start, end, offset, pos, rpos, rseg,
        #      segment_index, self.segment_list.nsegments,
        #      self.segment_list.segments[segment_index] )
        overlap = 1

        assert overlap > 0, \
            "sample %i-%i does not overlap with workspace: " \
            "offset=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
            (start, end, offset, rpos, rseg,
             segment_index, self.segment_list.nsegments,
             self.segment_list.segments[segment_index])

        return start, end, overlap

    def __dealloc__(self):
        if self.cdf != NULL:
            free(self.cdf)
            self.cdf = NULL


cdef class SegmentListSampler:
    '''return a new position within segment list.

    The probability of overlap is proportional
    to the workspace size.

    In order to avoid edge effects the position
    of the sample within a workspace segment is
    sampled again.
    '''

    cdef SegmentList segment_list
    cdef Position * cdf
    cdef Position total_size
    cdef int nsegments

    def __init__(self, SegmentList segment_list ):
        cdef Position i, totsize
        assert len(segment_list) > 0, "sampling from empty segment list"
        assert segment_list.isNormalized

        self.segment_list = segment_list
        self.nsegments = len(segment_list)
        self.cdf = <Position*>malloc( sizeof(Position) * self.nsegments )
        if not self.cdf:
            raise MemoryError(
                "out of memory when allocation %i bytes" %
                sizeof(sizeof(Position) * self.nsegments))
        self.total_size = 0
        for i from 0 <= i < len(segment_list):
            self.total_size += segment_length( segment_list.segments[i] )
            # -1, as searchsorted inserts on left
            self.cdf[i] = self.total_size -1

    cpdef sample( self, Position sample_length ):
        '''return a sampled segment of size sample_length.

        
        The length of the sampled segment might be less than *sample_length*
        if it extends beyond a workspace boundary.
        '''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef Position start, end
        cdef PositionDifference overlap
        cdef size_t segment_index
        cdef Segment chosen_segment

        cdef Position random_pos_in_workspace
        cdef long random_pos_in_segment, sampling_start
        # sample a position 
        # r = rand() / (RAND_MAX / total_size + 1)
        random_pos_in_workspace = numpy.random.randint(0, self.total_size)
        segment_index = searchsorted( self.cdf,
                                      self.nsegments,
                                      sizeof(Position),
                                      &random_pos_in_workspace,
                                      &cmpPosition,
                                      )

        assert 0 <= segment_index < self.nsegments, \
            "range error for %i: %i >= %i" % (random_pos_in_workspace, segment_index, self.nsegments)

        chosen_segment = self.segment_list.segments[segment_index]

        # resample start position within chosen segment
        # permitting for partial overlap.
        
        # start position allowed to be negative in order to avoid edge effects 
        # at the beginnig of a workspace 
        # if chosen_segment.start is < than sample_length
        sampling_start = <long>chosen_segment.start - sample_length + 1

        # reduce sampling space if previous segment is closer than sample_length
        # to avoid edge effects at workspace boundaries where the distance between
        # the adjacent workspaces is less than sample_length
        if segment_index > 0:
            sampling_start = lmax( self.segment_list.segments[segment_index-1].end,
                                   sampling_start)
        random_pos_in_segment = numpy.random.randint( \
            sampling_start,
            chosen_segment.end )

        # adjust position by sampled amount
        # adjust start position to negative coordinate
        # sampled length might be less sample_length
        start = lmax(0, random_pos_in_segment)
        end = random_pos_in_segment + sample_length

        # print "sample %i-%i, rand_pos_workspace=%i, rand_pos_segment=%i, index=%i/%i, segment=%s" %\
        #     (start, end, random_pos_in_workspace, random_pos_in_segment, 
        #      segment_index, self.segment_list.nsegments,
        #      self.segment_list.segments[segment_index] )
        overlap = range_overlap( chosen_segment.start, chosen_segment.end, start, end )

        assert overlap > 0, "sample %i-%i does not overlap with workspace: rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
            (start, end, random_pos_in_workspace, random_pos_in_segment,
             segment_index,
             self.segment_list.nsegments,
             self.segment_list.segments[segment_index])

        return start, end, overlap

    def __dealloc__(self):
        if self.cdf != NULL:
            free(self.cdf)
            self.cdf = NULL

cdef class HistogramSamplerSlow:

    cdef numpy.ndarray cdf
    cdef Position bucket_size

    def __init__(self,
                 numpy.ndarray histogram,
                 Position bucket_size):

        assert len(histogram) > 0, "sampling from empty histogram"

        self.cdf = numpy.cumsum(histogram, dtype = numpy.float)
        self.cdf /= self.cdf[-1]
        self.bucket_size = bucket_size

    cpdef Position sample( self ):
        '''return a new position within segment list.'''

        cdef Position base, ip
        cdef double r
        cdef Position bucket_size
        bucket_size = self.bucket_size
        # note: could be made quicker by
        # - creating a sample random doubles once, see numpy.random_integers?
        # - avoiding the binary search?
        r = numpy.random.random_sample( 1 )[0]
        ip = numpy.searchsorted( self.cdf, r )
        base = ip*bucket_size
        if bucket_size>1 : return base + numpy.random.randint(0, bucket_size)

        return base

cdef class HistogramSampler:

    cdef Position * cdf
    cdef Position bucket_size
    cdef Position nbuckets
    cdef Position total_size

    @cython.boundscheck(False)
    def __init__(self,
                 numpy.ndarray[DTYPE_INT_t, ndim=1] histogram,
                 Position bucket_size ):
        cdef Position i
        self.cdf = NULL

        self.nbuckets = len(histogram)
        assert self.nbuckets > 0, "sampling from empty histogram"

        self.cdf = <Position*>malloc( sizeof(Position) * self.nbuckets )
        if not self.cdf: raise MemoryError( "out of memory when allocation %i bytes" % sizeof( sizeof(Position) * self.nbuckets ))
        self.total_size = 0
        for i from 0 <= i < self.nbuckets:
            self.total_size += histogram[i]
            self.cdf[i] = self.total_size

        self.bucket_size = bucket_size

    cpdef Position sample( self ):
        '''return a new position within segment list.'''

        cdef Position base, index, r

        # 1 to avoid 0 length
        if self.total_size > 1:
            r = numpy.random.randint( 1, self.total_size )
        else:
            r = 1

        index = searchsorted( self.cdf,
                              self.nbuckets,
                              sizeof(Position),
                              &r,
                              &cmpPosition,
                              )
        base = index * self.bucket_size

        if self.bucket_size > 1 :
            return base + numpy.random.randint(0, self.bucket_size)

        return base

    def __dealloc__(self ):
        if self.cdf != NULL:
            free(self.cdf)
            self.cdf = NULL

cdef class Sampler:
    pass

cdef class SamplerAnnotator(Sampler):
    '''
    sample with the annotator method.

    Samples are created in a two step procedure. First, the length
    of a sample segment is chosen randomly from the empirical segment
    length distribution. Then, a random coordinate is chosen. The sampling
    stops when exactly the same number of nucleotides in the sampled segments 
    overlap the workspace as in the input.

    This method does not truncate segments to the workspace before
    or after sampling. The truncation can happen before and/or after
    the algorithm finishes.

    The sampling algorithm avoids checking for segment overlap until
    the number of nucleotides sampled reaches the expected count.
    Then, overlapping and adjacent segments are merged and the actual
    number of nucleotides overlapping the workspace computed.

    In case of over-shoot, a random segment and position within the current 
    sample is chosen and the overshoot is removed by trimming in a random
    direction.

    In case of under-shoot, the sampling process continues as usual. It 
    will converge if a segment of the appropriate length is sampled and
    placed without overlap or a sampled segment overlaps a previously
    sampled segment and due to the merging process increments its length
    by the correct number of bases.

    Conceptually this corresponds to a random process where the probability
    that a segment gains or loses a certain number of nucleotides is independent
    of the size of the segment (Poisson Process). An example of segments of this
    type are regions of conservation, indel-purified segments or inter-gap 
    segment lengths.

    Note that the length distribution of the sampled segments is not exactly 
    the same as the input segments. Also, the number of segments sampled will
    differ from the number of segments in the input set.

    This method is quick if the workspace is large compared to the segments that
    need to be placed. If it is small, a large number of samples will be rejected.

    *bucket_size* - bin size for empirical length distribution

    *nbuckets* - number of bins of the empirical length distribution

    The product of *bucket_size* and *nbuckets* needs to be larger than
    the largest segment in the input set.

    returns a SegmentList of sampled segments.
    '''


    cdef Position bucket_size
    cdef int nbuckets
    cdef int nunsuccessful_rounds

    def __init__( self,
                  bucket_size=1,
                  nbuckets=100000,
                  nunsuccessful_rounds=0):
        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.nunsuccessful_rounds = nunsuccessful_rounds

    def __reduce__(self):
        return (buildSamplerAnnotator, (self.bucket_size, 
                                        self.nbuckets, 
                                        self.nunsuccessful_rounds))

    cpdef SegmentList sample(self,
                             SegmentList segments,
                             SegmentList workspace):
        '''return a sampled list of segments.'''

        cdef PositionDifference remaining
        cdef PositionDifference true_remaining
        cdef PositionDifference overlap, ltotal, length

        cdef int nunsuccessful_rounds
        cdef int max_unsuccessful_rounds
        cdef Position start, end
        cdef Segment segment, sampled_segment
        cdef SegmentList sampled_segments
        cdef SegmentList unintersected_segments
        cdef SegmentList intersected_segments
        cdef SegmentList working_segments
        cdef SegmentList tmp_segments
        cdef SegmentListSampler sls, temp_sampler

        assert segments.isNormalized, "segment list is not normalized"
        assert workspace.isNormalized, "workspace is not normalized"

        unintersected_segments = SegmentList()
        intersected_segments = SegmentList()

        # collect all segments in workspace
        # This method does not truncate.
        working_segments = segments.clone()
        working_segments.filter(workspace)
        if len(working_segments) == 0:
            return intersected_segments

        # get nucleotides that need to be sampled
        # only count overlap within workspace
        tmp_segments = working_segments.clone()
        tmp_segments.intersect(workspace)
        ltotal = tmp_segments.sum()        

        # create space for sampled segments, add additional 10%
        # safety margin to avoid realloc calls
        sampled_segments = SegmentList(allocate=int(1.1 * len(segments)))

        # build a segment length histogram
        histogram, bucket_size = working_segments.getLengthDistribution(self.bucket_size,
                                                                        self.nbuckets)

        hs = HistogramSampler(histogram, bucket_size)

        # set up segment sampler
        sls = SegmentListSampler(workspace)

        remaining = ltotal
        true_remaining = remaining
        nunsuccessful_rounds = 0
        max_unsuccessful_rounds = 20

        while true_remaining > 0 and nunsuccessful_rounds < max_unsuccessful_rounds:

            ###################################################################
            # Sample a segment length from the histogram
            length = hs.sample()
            assert length > 0

            ###################################################################
            # If we can now potentially get required number of nucleotides, recompute
            # the required amount since we might need more due to overlap
            if remaining <= length:
                # add current sampled segments to list of segments
                unintersected_segments.extend(sampled_segments)
                # merge overlapping and adjacent segments
                unintersected_segments.merge(0)

                # reset sampled segments
                sampled_segments.clear()

                ############################################
                # compute overlap with workspace
                intersected_segments = unintersected_segments.clone()
                intersected_segments.intersect(workspace)

                assert intersected_segments.sum() <= unintersected_segments.sum()

                # recompute bases remaining to be covered
                remaining = ltotal - <PositionDifference>intersected_segments.sum()

                # check if we managed to increase coverage
                if true_remaining == remaining:
                    nunsuccessful_rounds += 1
                else:
                    true_remaining = remaining

            # deal with overshoot
            if true_remaining < 0:
                temp_sampler = SegmentListSampler( unintersected_segments )

                # sample a position from current list of sampled segments
                start, end, overlap = temp_sampler.sample(1)

                # causes a potential bias for overlapping segments:
                # in case of overlapping segments, the
                # segment chosen is always the first
                unintersected_segments.trim_ends(
                    start, 
                    -true_remaining,
                    numpy.random.randint(0, 2))

                # trimming might remove residues outside the workspace
                # hence - recompute true_remaining by going back to loop.
                true_remaining = 1
                continue

            # sample a position until we get a nonzero overlap
            start, end, overlap = sls.sample( length )

            sampled_segment = Segment( start, end )

            if true_remaining > 0:
                sampled_segments._add( sampled_segment )
                remaining -= overlap

        self.nunsuccessful_rounds = nunsuccessful_rounds

        # merge overlapping/adjacent
        intersected_segments = unintersected_segments.clone()
        intersected_segments.merge(0)

        # remove all segments outside workspace
        intersected_segments.filter(workspace)

        assert intersected_segments.sum() > 0
        return intersected_segments

def buildSamplerAnnotator(*args):
    '''unpickling - return a rebuild SamplerAnnotator object.'''
    return SamplerAnnotator(*args)


cdef class SamplerSegments(Sampler):
    '''
    sample *n* segments from segment length distribution.

    Samples are created in a two step procedure. First, the length
    of a sample segment is chosen randomly from the empirical segment
    length distribution. Then, a random coordinate is chosen. The sampling
    stops when exactly *n* segments have been placed.
    
    Sampled segments may overlap.

    The segment length distribution is derived from
    the argument *segments* supplied to the :meth:`sample`
    method.

    *n* is given by the number of segments in the observed
    data.

    *segments* does not need to be normalized. In fact,
    supplying overlapping segments is likely to be the 
    most realistic use case.

    *bucket_size* - bin size for empirical length distribution

    *nbuckets* - number of bins of the empirical length distribution

    The product of *bucket_size* and *nbuckets* needs to be larger than
    the largest segment in the input set.

    returns a SegmentList of sampled segments.
    '''

    cdef Position bucket_size
    cdef Position nbuckets

    def __init__( self, bucket_size = 1, nbuckets = 100000 ):

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets

    def __reduce__(self):
        return (buildSamplerSegments, (self.bucket_size, 
                                       self.nbuckets ))

    cpdef SegmentList sample( self,
                                           SegmentList segments,
                                           SegmentList workspace ):
        '''return a sampled list of segments.'''

        cdef Position length
        cdef SegmentListSampler sls
        cdef SegmentList sample

        assert workspace.isNormalized, "workspace is not normalized"

        sample = SegmentList(allocate=len(segments))
    
        # collect all segments in workspace
        working_segments = segments.clone()
        working_segments.filter(workspace)

        if len(working_segments) == 0:
            return sample

        # build length histogram
        histogram, bucket_size = working_segments.getLengthDistribution(
            self.bucket_size,
            self.nbuckets)
        hs = HistogramSampler(histogram, bucket_size)

        # create segment sampler
        sls = SegmentListSampler(workspace)

        for x in xrange(len(segments)):

            # Sample a segment length from the histogram
            length = hs.sample()
            assert length > 0

            # sample a position until we get a nonzero overlap
            start, end, overlap = sls.sample(length)

            sample._add(Segment(start, end))

        return sample

def buildSamplerSegments( *args ):
    '''unpickling - return a rebuild SamplerAnnotator object.'''
    return SamplerAnnotator( *args )

########################################################
########################################################
########################################################
cdef class SamplerBruteForce(Sampler):
    '''
    sample segments from length distribution until
    workspace is covered by at least the same number of nucleotides
    as the input.

    This sampler applies a brute force algorithm - segments 
    are added only if the do not overlap with a previously
    sampled segment.

    The segment length distribution is derived from
    the argument *segments* supplied to the :meth:`sample`
    method.

    *bucket_size* - bin size for empirical length distribution

    *nbuckets* - number of bins of the empirical length distribution

    The product of *bucket_size* and *nbuckets* needs to be larger than
    the largest segment in the input set.

    returns a SegmentList of sampled segments.

    '''

    cdef Position bucket_size
    cdef Position nbuckets
    cdef int ntries_inner
    cdef int ntries_outer

    def __init__( self, 
                  bucket_size=1, 
                  nbuckets=100000,
                  ntries_inner=100,
                  ntries_outer=10):

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.ntries_inner = ntries_inner
        self.ntries_outer = ntries_outer

    def __reduce__(self):
        return (buildSamplerBruteForce, (self.bucket_size, 
                                         self.nbuckets, 
                                         self.ntries_inner,
                                         self.ntries_outer))

    cpdef SegmentList sample(self,
                             SegmentList segments,
                             SegmentList workspace):
        '''return a sampled list of segments.'''

        cdef Position start, end, length
        cdef SegmentListSampler sls
        cdef SegmentList sample
    
        assert workspace.isNormalized, "workspace is not normalized"

        sample = SegmentList(allocate=len(segments))

        # collect all segments in workspace
        working_segments = segments.clone()
        working_segments.filter(workspace)

        if len(working_segments) == 0:
            return sample

        # build length histogram
        histogram, bucket_size = working_segments.getLengthDistribution(
            self.bucket_size,
            self.nbuckets)
        hs = HistogramSampler(histogram, bucket_size)

        # create segment sampler
        sls = SegmentListSampler(workspace)

        cdef int ntries_outer = self.ntries_outer
        cdef int ntries_inner
        cdef PositionDifference remaining, overlap

        while ntries_outer > 0:

            sample.clear()

            remaining = segments.sum()
            ntries_inner = self.ntries_inner
            # print "starting inner tries"

            while remaining > 0 and ntries_inner > 0:

                # Sample a segment length from the histogram
                length = hs.sample()
                assert length > 0

                # sample a position until we get a nonzero overlap
                start, end, overlap = sls.sample(length)

                # print str(sample), start, end, overlap, remaining

                if overlap > remaining:
                    ntries_inner -= 1
                    continue

                if sample.overlapWithRange(start, end): 
                    ntries_inner -= 1
                    continue

                sample._add(Segment(start, end))
                # required in order to sort samples
                sample.normalize()

                # print "adding", start, end, str(sample)

                ntries_inner = self.ntries_inner

                remaining -= overlap

            if ntries_inner > 0:
                break

            ntries_outer -= 1
            
        if ntries_outer == 0:
            raise ValueError("sampling did not converge: %s" % str(sample))

        return sample

def buildSamplerBruteForce(*args):
    '''unpickling - return a rebuild SamplerAnnotator object.'''
    return SamplerBruteForce(*args)


cdef class SamplerUniform(Sampler):
    '''
    For debugging purposes.

    Every *increment* residues within the :term:`workspace`
    a new segment is generated. The segment length is
    taken from the segment size distribution of the input set.

    Sampled segments extend from the position selected alternately into 
    into forward or reverse direction.

    *bucket_size* - bin size for empirical length distribution

    *nbuckets* - number of bins of the empirical length distribution

    The product of *bucket_size* and *nbuckets* needs to be larger than
    the largest segment in the input set.
    '''

    cdef Position bucket_size
    cdef Position nbuckets
    cdef Position increment
    cdef Position current_workspace
    cdef Position current_position
    cdef int current_orientation

    def __init__( self, 
                  increment,
                  bucket_size = 1,
                  nbuckets = 100000,
                  current_orientation = 0,
                  current_workspace = 0,
                  current_position = 0):
        self.current_orientation = current_orientation
        self.current_workspace = current_workspace
        self.current_position = current_position

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.increment = increment

    def __reduce__(self):
        return (buildSamplerUniform, (self.increment,
                                      self.bucket_size, 
                                      self.nbuckets, 
                                      self.current_orientation,
                                      self.current_workspace,
                                      self.current_position) )

    cpdef SegmentList sample( self,
                                           SegmentList segments,
                                           SegmentList workspace ):
        '''return a sampled list of segments.'''
        assert workspace.isNormalized, "workspace is not normalized"

        cdef SegmentList sample
        cdef SegmentListSampler sls
        cdef Position increment = self.increment
        cdef Position start, end, x, length, i, added, nsegments, nworkspaces
        cdef Position current_offset, current_workspace

        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )
        working_segments.filter( workspace )

        # allocate sample large enough
        sample = SegmentList( allocate = (len(workspace) / increment + workspace.sum() ) )

        if len(working_segments) == 0:
            return sample

        # build length histogram
        histogram, bucket_size = working_segments.getLengthDistribution( self.bucket_size,
                                                                         self.nbuckets )

        hs = HistogramSampler( histogram, bucket_size )

        sample = SegmentList( allocate = increment )

        added = 0

        nsegments = len(working_segments)
        nworkspaces = len(workspace)

        x = self.current_position 
        current_workspace = self.current_workspace
        start, end = workspace[current_workspace]        

        while added < nsegments:

            while x > end:
                x -= end
                current_workspace = (current_workspace + 1) % nworkspaces
                start, end = workspace[current_workspace]        
                x += start            

            # sample a segment length from the histogram
            length = hs.sample()
            assert length > 0
            if self.current_orientation:
                sample._add( Segment( x, x + length ) )
                self.current_orientation = 0
            else:
                sample._add( Segment( x-length, x) )
                self.current_orientation = 1
            x += increment
            added += 1
                             
        self.current_position = x
        self.current_workspace = current_workspace

        return sample

def buildSamplerUniform( *args ):
    '''unpickling - return a rebuilt SamplerUniform object.'''
    return SamplerUniform( *args )

########################################################
########################################################
########################################################
cdef class SamplerShift(Sampler):
    '''
    Sample segments by shifting them by a random amount within *radius*
    in a random direction.

    This sampler can be used to investigate overlap between
    two types of genomic annotations that tend to co-occur.
    An example are two transcription factors that tend to 
    occur in the promotor of the same genes. 

    *radius* determines the size of the region that is accessible 
    for randomly shifting a segment. It is expressed as a fraction 
    of the size of a segment.

    In the sampling procedure, the start of the segment is shifted
    within radius. It is then wrapped around any discontinuities in the workspace. If
    the segment extends beyond the radius, the remaining nucleotides will
    be wrapped around.
    '''

    cdef double radius
    cdef int extension 

    def __init__( self, 
                  radius = 2,
                  extension = 0):

        self.radius = radius
        self.extension = extension

    def __reduce__(self):
        return (buildSamplerShift, (self.radius, 
                                    self.extension)), 
    
    cpdef SegmentList sample(self,
                             SegmentList segments,
                             SegmentList workspace):
        '''create a random sample of segments.

        *segments* - a list of segments
        
        *workspace* - the workspace to use

        returns: a list of sampled segments.
        '''

        assert workspace.isNormalized, "workspace is not normalized"

        cdef Segment segment
        cdef SegmentList sample, working_segments
        cdef Position length, direction, extended_length
        cdef Position x, midpoint
        cdef PositionDifference shift, start, end, ws_start, ws_end, shift_area
        cdef PositionDifference remainder, remainder_left, remainder_right
        cdef double half_radius = self.radius / 2
        cdef PositionDifference half_extension = self.extension // 2
        cdef Segment * _working_segments
        cdef SegmentList ws

        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )
        working_segments.filter( workspace )
        _working_segments = working_segments.segments

        # allocate sample large enough
        sample = SegmentList( len(working_segments) * 2 )

        if len(working_segments) == 0:
            return sample
        
        for x from 0 <= x < len(working_segments):
            segment = _working_segments[x]
            length = segment.end - segment.start
            midpoint = segment.start + length // 2
            if self.extension:
                shift_area = half_extension
            else:
                shift_area = <Position>floor(length * half_radius)

            # get workspace around segment
            ws_start = lmax(0, midpoint - shift_area)
            ws_end = lmax(0, midpoint + shift_area) 
            ws = workspace.getOverlappingSegmentsWithRange( ws_start, ws_end )
            ws.truncate( Segment( ws_start, ws_end ) )

            start = ws.getRandomPosition()
            if numpy.random.randint( 0,2 ):
                end = start + length
            else:
                end = start
                start = end - length

            # ws_start is now the intersection of workspace and sampling space
            ws_start, ws_end = ws.min(), ws.max()

            remainder = length
            if start < ws_start:
                # start might be further than length outside of range
                remainder = lmin(ws_start - start, length)
                sample.extend( ws.getFilledSegmentsFromStart( start, length - remainder ) )
                sample.extend( ws.getFilledSegmentsFromEnd( ws_end, remainder ) )
            elif end > ws_end:
                # end might be further than length outside of range
                remainder = lmin(end - ws_end, length)
                sample.extend( ws.getFilledSegmentsFromEnd( end, length - remainder ) )
                sample.extend( ws.getFilledSegmentsFromStart( ws_start, remainder ) )
            else:
                sample.extend( ws.getFilledSegmentsFromStart( start, length ) )

        sample.normalize()
        return sample

def buildSamplerShift( *args ):
    '''unpickling - return a rebuilt SamplerShift object.'''
    return SamplerShift( *args )

########################################################
########################################################
########################################################
cdef class SamplerLocalPermutation(Sampler):
    '''
    Sample segments by local permutation.

    This sampler creates a randomized copy of the input
    set of segmeents by local permutation. The sampler proceeds locally,
    that is for each continuous part of the workspace separately.

    First, :term:`segments` of the input set overlapping with the current section
    of the workspace are collected. Then, the order of the segments is randomly
    permuted and randomly sized gaps inserted. Finally, segments are entered 
    into the workspace from a randomly chosen point. If a segment extends beyond
    the workspace boundary, it is wrapped around to the start of the
    workspace.

    If a segment of the input set is extending beyond a workspace boundary, the
    full sized segment is used for the permutation, but the workspace segment
    is enlarged by the same amount::

           |-----workspace segment--------------------|
      1111111      222         333      44444 55555
      |----------workspace size for permutation-------|
      33 4444    1111111     55555         222        3   - Random sample 1
        222     44444 1111111    333   55555              - Random sample 2


    '''

    def __init__( self ):
        pass

    cpdef SegmentList sample( self,
                                           SegmentList segments,
                                           SegmentList workspace ):
        '''create a random sample of segments.

        .. note::
            Still needs to optimized for speed.
            The algorithm is not fully implemented yet.

        *segments* - a list of segments
        
        *workspace* - the workspace to use

        returns: a list of sampled segments.
        '''

        assert workspace.isNormalized, "workspace is not normalized"

        cdef Segment segment
        cdef SegmentList working_segments
        cdef PositionDifference work_start, work_end, start, end, last
        cdef PositionDifference total_length, free_length
        cdef Segment* _working_segments
        cdef PositionDifference shift
        cdef SegmentList sample = SegmentList()
        
        for work_start, work_end in workspace:

            # collect all segments in this workspace segment
            working_segments = segments.getOverlappingSegments( Segment( work_start, work_end) )
            if len(working_segments) == 0: continue
            _working_segments = working_segments.segments
            
            lengths = [ x[1] - x[0] for x in working_segments ]
            total_length = sum( lengths )
            
            # modify workspace to deal with overhanging segments
            # extend workspace
            work_start = lmin( working_segments.min(), work_start )
            work_end = lmax( working_segments.max(), work_end )
            free_length = work_end - work_start - total_length
            
            # 1. permutate order of segments
            random.shuffle( lengths )
        
            # 2. determine size of space between samples
            # sample points in free area. The distance
            # between the points is the space left.
            points = []
            for x in range(len(lengths) ):
                points.append( random.randint( 0, free_length ) )
            points.sort()

            # cycle shift to avoid edge effects
            shift = random.randint( 0, free_length )

            # 3. move segments to appropriate place
            start = work_start + shift
            last = 0
            for x in range(len(lengths)):
                start += points[x] - last
                # wrap around if beyond end of workspace
                if start > work_end:
                    start = work_start + start - work_end
                end = start + lengths[x]
                if end < work_end:
                    sample._add( Segment(start, end ) )
                    start += lengths[x]
                else:
                    # wrap around segment
                    sample._add( Segment(start, work_end ) )
                    end = work_start + end - work_end
                    sample._add( Segment( work_start, end ) )

                start = end
                last = points[x]
                
            assert start + (points[-1] - last) <= work_end, "start=%i, points[-1]=%i, work_end=%i" % \
                (start, points[-1] -last, work_end)
        
            sample.normalize()
        return sample

########################################################
########################################################
########################################################
cdef class SamplerGlobalPermutation(Sampler):
    '''
    Sample segments by global permutation.

    This sampler creates a randomized copy of the input
    set of segmeents by permutation. The sampler uses all segments
    within the workspace.

    First, :term:`segments` of the input set overlapping with the workspace
    are collected. Then, the order of the segments is randomly
    permuted and randomly sized gaps inserted. Finally, segments are entered 
    into the workspace from a randomly chosen point. If a segment extends beyond
    the end of workspace segment it is wrapped around to the start of the next
    workspace segment. 

    If a segment of the input set is extending beyond a workspace boundary, the
    full sized segment is used for the permutation, but the workspace segment
    is enlarged by the same amount::

           |--workspace segment--|      |--workspace segment--|
      1111111      222         333      44444 55555
      |--------------------------|      |---------------------|  - Workspace used for permutation
      33 4444    1111111       55       555         222       3  - Random sample 1
        222     44444 1111111              333   55555           - Random sample 2

    '''

    def __init__( self ):
        pass

    cpdef SegmentList sample( self,
                                           SegmentList segments,
                                           SegmentList workspace ):
        '''create a random sample of segments.

        .. note::
            Still needs to optimized for speed

        *segments* - a list of segments
        
        *workspace* - the workspace to use

        returns: a list of sampled segments.
        '''

        assert workspace.isNormalized, "workspace is not normalized"

        cdef Segment segment
        cdef SegmentList working_segments, working_workspace
        cdef PositionDifference work_start, work_end, start, end, last_end
        cdef PositionDifference total_length, free_length
        
        cdef SegmentList sample = SegmentList()

        # collect all segments in workspace
        # This method does not truncate.
        working_segments = SegmentList( clone = segments )
        working_segments.filter( workspace )

        if len(working_segments) == 0: return sample
        
        # create a copy of the workspace
        working_workspace = SegmentList( clone = workspace )
        
        # extend workspace segments with segments
        working_workspace.extend( working_segments )
        working_workspace.merge(0)

        # get lengths of segments inside and extending from workspace
        lengths = [ x[1] - x[0] for x in working_segments ]
        total_length = sum( lengths )

        # compute workspace size free of segments
        free_length = working_workspace.sum()
        free_length -= total_length

        # 1. permutate order of segments
        random.shuffle( lengths )

        # 2. determine size of space between samples
        # sample points in free area. The distance
        # between the points is the space left.
        points = []
        for x in range(len(lengths) ):
            points.append( random.randint( 0, free_length ) )
        points.sort()

        # 3. cycle shift to avoid edge effects
        shift = random.randint( 0, free_length )

        # find starting position 
        count = 0
        idx = 0
        for start,end in working_workspace:
            count += end - start
            if count > shift: 
                count -= end - start
                break
            idx += 1

        # idx is now the segment to start inserting in
        # add missing residues
        start += shift - count

        # 4. move segments to appropriate place
        last_end = 0
        global_work_start, global_work_end = working_workspace.min(), working_workspace.max()
        work_start, work_end = working_workspace[idx]

        #print "starting adding", start
        max_idx = len(working_workspace)
        last_ungapped_position = 0
        for segment_idx in range( len(lengths) ):
            
            # place a gap
            gap_length = points[segment_idx] - last_ungapped_position
            assert gap_length >= 0
            #print "gap_length", gap_length
            # print "gap"
            while gap_length > 0:
                increment = lmin( work_end - start, gap_length )
                # print idx, start, start + increment, gap_length
                start += increment
                if start == work_end:
                    idx += 1
                    if idx >= max_idx: idx = 0
                    work_start, work_end = working_workspace[idx]
                    start = work_start
                gap_length -= increment
                last_ungapped_position += increment

            # place a segment of a certain length
            # the segment is split into multiple components at
            # workspace gaps
            length = lengths[segment_idx]
            #print "segment", length

            while length > 0:
                increment = lmin( work_end - start, length )
                #print idx, 'start=',start, 'end=', start + increment, 'len=',length
                if increment > 0:
                    sample._add( Segment( start, start + increment ) )
                start += increment
                if start == work_end:
                    idx += 1
                    if idx >= max_idx: idx = 0
                    work_start, work_end = working_workspace[idx]
                    start = work_start
                length -= increment

        sample.normalize()

        return sample

########################################################
########################################################
########################################################
cdef class SamplerDummy(Sampler):

    def __init__( self ):
        '''
        returns a copy of the observed data.
        '''

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''

        cdef SegmentList sample
        sample = SegmentList( clone = segments )
        return sample

############################################################
############################################################
############################################################
## Counters
############################################################
cdef class Counter:
    '''base class for objects that compute counts
    between two segment lists.
    '''

cdef class CounterNucleotideOverlap(Counter):
    """return number of nucleotides overlapping between segments and
    annotations."""

    name = "nucleotide-overlap"

    def __call__(self, SegmentList segments,
                 SegmentList annotations,
                 SegmentList workspace=None):
        return annotations.overlapWithSegments(segments)

cdef class CounterNucleotideDensity(Counter):
    name = "nucleotide-density"

    def __call__(self, SegmentList segments,
                 SegmentList annotations,
                 SegmentList workspace ):
        '''return number of nucleotides overlapping between segments and annotations.
        divided by the size of the workspace
        '''
        cdef Position l
        l = len(workspace)
        if l == 0:
            return 0
        return float(annotations.overlapWithSegments(segments)) / l

cdef class CounterSegmentOverlap(Counter):
    name = "segment-overlap"

    def __call__(self, segments, annotations, workspace=None):
        '''return number of segments overlapping with annotations.'''
        return segments.intersectionWithSegments(annotations)

cdef class CounterSegmentMidpointOverlap(Counter):
    name = "segment-midoverlap"

    def __call__(self, segments, annotations, workspace=None ):
        '''return number of segments overlapping with annotations.'''
        return segments.intersectionWithSegments(annotations,
                                                 mode="midpoint")

cdef class CounterAnnotationOverlap(Counter):
    name = "annotation-overlap"

    def __call__(self, segments, annotations, workspace=None):
        '''return number of annotations overlapping with segments.'''
        return annotations.intersectionWithSegments(segments)

cdef class CounterAnnotationMidpointOverlap(Counter):
    name = "annotation-midoverlap"

    def __call__(self, segments, annotations, workspace=None):
        '''return number of annotations overlapping with segments.'''
        return annotations.intersectionWithSegments(
            segments,
            mode="midpoint")

############################################################
############################################################
############################################################
## Annotator results
############################################################
@cython.boundscheck(False)
cpdef getNPTwoSidedPValue(ar, val):
    '''return pvalue for *val* within sorted array *ar*
    '''
    idx = numpy.searchsorted(ar, val)
    l = len(ar)
    min_pval = 1.0 / l
    if idx == l:
        pval = min_pval
    elif idx > l / 2:
        # over-representation
        while idx > 0 and ar[idx] == val: idx -= 1
        pval = 1.0 - float(idx) / l
    else:
        # under-representation
        while idx < l and ar[idx] == val: idx += 1
        pval = float(idx) / l

    return dmax(min_pval, pval)

@cython.boundscheck(False)
def getNPTwoSidedPValueFast(numpy.ndarray[DTYPE_FLOAT_t, ndim=1] ar,
                            double val,
                            double mean):
    '''return a two-sided pvalue.

    Fast if val is small or large.
    '''
    cdef Position l, x
    cdef double min_pval, pval

    l = len(ar)
    min_pval = 1.0 / l

    if val > mean:
        x = 0
        while x < l and ar[x] < val: x += 1
        pval = float(x) / l
    else:
        x = l - 1
        while x >= 0 and ar[x] > val: x -= 1
        pval = 1.0 - float(x) / l

    return dmax(min_pval, pval)

############################################################
############################################################
############################################################
## Annotator results
############################################################
ctypedef struct EnrichmentStatistics:
    double expected
    double observed
    double stddev
    double lower95
    double upper95
    double fold
    Position nsamples
    double * samples
    int * sorted2sample
    int * sample2sorted
    double pvalue
    double qvalue 

cdef double getTwoSidedPValue( EnrichmentStatistics * stats, 
                               double val ):
    '''return pvalue for *val* within sorted array *ar*
    '''
    cdef long idx, l
    cdef double min_pval, pval, rval

    # find index of value correspending to observed value in samples
    idx = searchargsorted( stats.samples,
                           stats.sorted2sample,
                           stats.nsamples,
                           sizeof(double),
                           &val,
                           &cmpDouble,
                        )
    
    l = stats.nsamples
    min_pval = 1.0 / l

    if idx == l:
        idx = 1
    elif val > stats.expected :
        # over-representation
        while idx > 0 and stats.samples[stats.sorted2sample[idx]] == val: idx -= 1
        idx = l - (idx+1)
    else:
        # under-representation
        while idx < l and stats.samples[stats.sorted2sample[idx]] == val: 
            idx += 1
        # no -1 because of 0-based indices

    pval = float(idx) / l
    
    return dmax( min_pval, pval)

cdef void compressSampleIndex( EnrichmentStatistics * stats ):
    '''compress indices in stats.'''

    cdef int x, refidx, observed_idx
    cdef double lastval, thisval
    cdef Position l
    l = stats.nsamples

    # normalize - equal values will get the same index
    # before
    # samples        0 0 0 1 1 1
    # sample2sorted  0 1 2 3 4 5
    # sorted2sample  0 1 2 3 4 5
    # after
    # samples        0 0 0 1 1 1
    # sample2sorted  2 2 2 3 3 3
    # sorted2sample  0 1 2 3 4 5
    #
    # for values < observed: last index
    # for values > observed: first index

    # locate midpoint differentiating over and under-represneted
    # observed_idx = index of element with first sample > observed
    observed_idx = 0
    while observed_idx < l and stats.samples[stats.sorted2sample[observed_idx]] <= stats.observed:
        observed_idx += 1

    # print "obs_idx=", observed_idx, "observed=", observed
    x = observed_idx - 1
    lastval = stats.samples[stats.sorted2sample[x]]
    refidx = x

    while x >= 0:
        # print "l", x
        thisval = stats.samples[stats.sorted2sample[x]]

        if thisval != lastval:
            lastval = thisval
            refidx = x
            
        stats.sample2sorted[stats.sorted2sample[x]] = refidx
        x -= 1 

    x = observed_idx
    lastval = stats.samples[stats.sorted2sample[x]]
    refidx = x

    while x < l:
        # print "r", x
        thisval = stats.samples[stats.sorted2sample[x]]
        if thisval != lastval:
            lastval = thisval
            refidx = x
            
        stats.sample2sorted[stats.sorted2sample[x]] = refidx
        x += 1

cdef EnrichmentStatistics * makeEnrichmentStatistics(observed,
                                                     samples,
                                                     reference,
                                                     pseudo_count):

    cdef EnrichmentStatistics * stats 
    cdef Position offset, i, l

    l = len(samples)
    if l < 1:
        return NULL

    stats = <EnrichmentStatistics*>malloc(sizeof(EnrichmentStatistics))
    if not stats:
        raise MemoryError("out of memory when allocation %i bytes" % sizeof(EnrichmentStatistics) )

    stats.samples = <double*>calloc( l, sizeof(double))
    stats.sorted2sample = <int*>calloc( l, sizeof(int))
    stats.sample2sorted = <int*>calloc( l, sizeof(int))

    stats.observed = observed
    stats.nsamples = l 
    for i from 0 <= i < l: stats.samples[i] = float(samples[i])

    # create index of sorted values
    r = numpy.argsort( samples )

    # save map: sample_id to sort order
    for i from 0 <= i < l: 
        stats.sample2sorted[r[i]] = i
        stats.sorted2sample[i] = r[i]

    # compressSampleIndex( stats )
    # print "after"
    # for i from 0 <= i < l: 
    #     print i, stats.samples[i], stats.sorted2sample[i], stats.sample2sorted[i]
    stats.expected = numpy.mean(samples)

    if reference != None:
        # test against reference
        # move expected by fold change in the reference set
        stats.expected *= reference.fold

    # optionally add pseudo_counts
    if stats.expected != 0:
        stats.fold = (stats.observed + pseudo_count) / (stats.expected + pseudo_count)
    else:
        stats.fold = 1.0

    stats.stddev = numpy.std(samples)

    # compute 95% confidence intervals 
    # The confidence interval are the values that lie
    # at the 5% and 95% percentile of the samples.
    offset = int(0.05 * l)
    
    if offset > 0: 
        stats.lower95 = stats.samples[stats.sorted2sample[ lmin( offset, l-1) ]]
        stats.upper95 = stats.samples[stats.sorted2sample[ lmax( l-offset, 0) ]]
    else:
        stats.lower95 = stats.samples[stats.sorted2sample[ 0 ]]
        stats.upper95 = stats.samples[stats.sorted2sample[ l-1 ]]
    
    # adjust confidence intervals for reference fold change
    # NB: I am not sure that this is proper
    if reference == None:
        stats.pvalue = getTwoSidedPValue( stats, stats.observed )
    else:
        # compute adjusted pvalue. Conceptually, move the distribution to center around 
        # expected instead of the mean. Instead of moving the distribution, simply changing 
        # the observed value is equivalent
        if reference.fold > 0:
            stats.pvalue = getTwoSidedPValue( stats, stats.observed / reference.fold )
        else:
            raise ValueError( "0 fold change not applicable" )

        # this is analogous to shifting the excepted value.
        # Is this correct?
        stats.lower95 *= reference.fold
        stats.upper95 *= reference.fold

    stats.qvalue = 1.0

    return stats

############################################################
############################################################
############################################################
## Annotator results
############################################################
cdef class AnnotatorResult(object):
    '''container for annotator results.'''

    cdef str format_observed
    format_expected = "%6.4f"
    format_fold = "%6.4f"
    format_pvalue = "%6.4e"
    format_counts = "%i"
    format_density = "%6.4e"

    headers = ["track", 
               "annotation",
               "observed",
               "expected",
               "CI95low", 
               "CI95high",
               "stddev",
               "fold",
               "l2fold",
               "pvalue",
               "qvalue",
               ]

    cdef:
        EnrichmentStatistics * stats
        str track, annotation, counter

    def __init__( self,
                  track,
                  annotation,
                  counter,
                  observed,
                  samples,
                  reference = None,
                  pseudo_count = 1.0 ):
        self.track = track
        self.annotation = annotation
        self.counter = counter
        self.stats = makeEnrichmentStatistics(observed, 
                                              samples,
                                              reference,
                                              pseudo_count)

        self.format_observed = "%i"

    def __str__(self):

        # if self.stats.nsamples < 10**6:
        #     format_pvalue = "%7.6f"
        # else:
        #     format_pvalue = "%7.6e"

        if self.stats.fold > 0:
            logfold = self.format_fold % math.log(self.stats.fold, 2)
        else:
            logfold = "-inf"

        return "\t".join( (self.track,
                           self.annotation,
                           self.format_observed % self.stats.observed,
                           self.format_expected % self.stats.expected,
                           self.format_expected % self.stats.lower95,
                           self.format_expected % self.stats.upper95,
                           self.format_expected % self.stats.stddev,
                           self.format_fold % self.stats.fold,
                           logfold,
                           self.format_pvalue % self.stats.pvalue,
                           self.format_pvalue % self.stats.qvalue,
                           ) )

    def __dealloc__(self):
        if self.stats != NULL:
            free(self.stats.samples)
            free(self.stats.sorted2sample)
            free(self.stats.sample2sorted)
            free(self.stats)

    property track:
        def __get__(self): return self.track

    property annotation:
        def __get__(self): return self.annotation

    property counter:
        def __get__(self): return self.counter

    property observed:
        def __get__(self): return self.stats.observed

    property expected:
        def __get__(self): return self.stats.expected

    property fold:
        def __get__(self): return self.stats.fold

    property stddev:
        def __get__(self): return self.stats.stddev

    property pvalue:
        def __get__(self): return self.stats.pvalue
        def __set__(self, val): self.stats.pvalue = val

    property qvalue:
        def __get__(self): return self.stats.qvalue
        def __set__(self, val): self.stats.qvalue = val
 
    property nsamples:
        def __get__(self): return self.stats.nsamples

    property samples:
        def __get__(self): 
            cdef Position x
            r = numpy.zeros( self.nsamples, dtype = numpy.float )
            for x from 0 <= x < self.stats.nsamples:
                r[x] = self.stats.samples[x]
            return r

    property format_observed:
        def __set__(self,f): self.format_observed = f

    def isSampleSignificantAtPvalue( self, sample_id, double pvalue ):
        return isSampleSignificantAtPvalue( self.stats, sample_id, pvalue )

    def getSample( self, sample_id ):
        return self.stats.samples[sample_id]

    def getEmpiricalPValue( self, value ):
        return getTwoSidedPValue( self.stats, value )

cdef class AnnotatorResultExtended(AnnotatorResult):
    '''container for annotator results.'''

    format_density = "%6.4e"

    headers = ["track", 
               "annotation",
               "observed",
               "expected",
               "CI95low", 
               "CI95high",
               "stddev",
               "fold",
               "l2fold",
               "pvalue",
               "qvalue",
               "track_nsegments",
               "track_size",
               "track_density",
               "annotation_nsegments",
               "annotation_size",
               "annotation_density",
               "overlap_nsegments",
               "overlap_size",
               "overlap_density",
               "percent_overlap_nsegments_track",
               "percent_overlap_size_track",
               "percent_overlap_nsegments_annotation",
               "percent_overlap_size_annotation",
               ]

    cdef:
        Position track_nsegments
        Position track_size
        Position annotation_nsegments
        Position annotation_size
        Position overlap_nsegments
        Position overlap_size
        Position workspace_size

    def __init__(self,
                 track,
                 annotation,
                 counter,
                 observed,
                 samples,
                 track_segments,
                 annotation_segments,
                 workspace,
                 reference=None,
                 pseudo_count=1.0):

        AnnotatorResult.__init__(
            self, track, annotation, counter, observed, samples, 
            reference = reference, 
            pseudo_count = pseudo_count)

        self.track_nsegments = track_segments.counts()
        self.track_size = track_segments.sum()

        self.annotation_nsegments = annotation_segments.counts()
        self.annotation_size = annotation_segments.sum()

        overlap = track_segments.clone()
        try:
            overlap.intersect(annotation_segments)
        except TypeError:
            # TODO: intersection between SegmentList with PositionList
            # needs still to be implemented.
            pass

        self.overlap_nsegments = overlap.counts()
        self.overlap_size = overlap.sum()

        self.workspace_size = workspace.sum()

    def __str__(self):

        # if self.stats.nsamples < 10**6:
        #     format_pvalue = "%7.6f"
        # else:
        #     format_pvalue = "%7.6e"

        if self.stats.fold > 0:
            logfold = self.format_fold % math.log( self.stats.fold, 2 )
        else:
            logfold = "-inf"

        def _toFold( a, b ):
            if b > 0: return self.format_fold % (100.0 * float(a) / b )
            else: return "na"

        def _toDensity( a, b ):
            if b > 0: return self.format_density % (100.0 * float(a) / b )
            else: return "na"

        return "\t".join( (self.track,
                           self.annotation,
                           self.format_observed % self.stats.observed,
                           self.format_expected % self.stats.expected,
                           self.format_expected % self.stats.lower95,
                           self.format_expected % self.stats.upper95,
                           self.format_expected % self.stats.stddev,
                           self.format_fold % self.stats.fold,
                           logfold,
                           self.format_pvalue % self.stats.pvalue,
                           self.format_pvalue % self.stats.qvalue,
                           self.format_counts % self.track_nsegments,
                           self.format_counts % self.track_size,
                           _toDensity( self.track_size, self.workspace_size),
                           self.format_counts % self.annotation_nsegments,
                           self.format_counts % self.annotation_size,
                           _toDensity( self.annotation_size, self.workspace_size),
                           self.format_counts % self.overlap_nsegments,
                           self.format_counts % self.overlap_size,
                           _toDensity( self.overlap_size, self.workspace_size),
                           _toFold( self.overlap_nsegments, self.track_nsegments ),
                           _toFold( self.overlap_size, self.track_size ),
                           _toFold( self.overlap_nsegments, self.annotation_nsegments ),
                           _toFold( self.overlap_size, self.annotation_size ),
                           ) )

############################################################
############################################################
############################################################
def getNormedPValue( value, r ):
    '''return pvalue assuming that samples are normal distributed.'''
    absval = abs(value - r.expected)
    if r.stddev == 0: 
        # dummy pvalue - if no overlap in expected, pvalue is set to 1.
        pvalue = 1.0
    else:
        if HAS_SCIPY:
            pvalue = 1.0 - scipy.stats.norm.cdf( absval, 0, r.stddev )
        else:
            raise ImportError( "scipy required" )
    return pvalue

############################################################
############################################################
############################################################
def getEmpiricalPValue( value, r ):
    return r.getEmpiricalPValue( value )

############################################################
############################################################
############################################################
def updatePValues( annotator_results, method = "empirical" ):
    '''update pvalues.

    empirical
        report pvalues from simulations. Minimum pvalue is 
        1/nsamples
    norm
        fit Gaussian to simulated values and compute pvalue from
        the distribution
    '''

    if method == "norm":
        methodf = getNormedPValue
    elif method == "empirical":
        methodf = getEmpiricalPValue
    else:
        raise ValueError( "unknown method '%s'" % method )

    for r in annotator_results:
        r.pvalue = methodf( r.observed, r )

############################################################
############################################################
############################################################
def getQValues( pvalues, method = "storey", **kwargs ):
    '''return a list of qvalues for a list of pvalues.'''
    
    if method == "storey":
        try:
            fdr = gat.Stats.computeQValues( pvalues, 
                                            vlambda = kwargs.get( "vlambda", numpy.arange( 0,0.95,0.05) ),
                                            pi0_method = kwargs.get( "pi0_method", "smoother") )
        except ValueError, msg:
            E.warn( "qvalue computation failed: %s" % msg )
            return [1.0] * len(pvalues)
        
        return fdr.qvalues
    else:
        return gat.Stats.adjustPValues( pvalues, method = method )

############################################################
############################################################
############################################################
def updateQValues( annotator_results, method = "storey", **kwargs ):
    '''update qvalues in annotator results

    storey
        qvalue from the method by Storey et al.
    '''

    pvalues = [ r.pvalue for r in annotator_results ]

    for r, qvalue in zip( annotator_results, getQValues( pvalues, method, **kwargs )):
        r.qvalue = qvalue

############################################################
############################################################
############################################################
## Workspace generators
############################################################
class UnconditionalWorkspace:
    '''compute a conditional workspace.

    the default workspace is not conditional.'''

    is_conditional = False

    def __call__(self, segments, annotations, workspace):
        return segments, annotations, workspace

    def filter( self, segments, annotations, workspace ):
        '''restrict annotations and segments to workspace.

        This speeds up computations.
        
        returns filtered segments, filtered annotations, workspace
        '''

        if annotations:
            temp_annotations = annotations.clone()
            temp_annotations.filter( workspace )
        else:
            temp_annotations = None
        
        if segments:
            temp_segments = segments.clone()
            temp_segments.filter(workspace)
        else:
            temp_segments = None

        return temp_segments, temp_annotations, workspace

class ConditionalWorkspaceCooccurance( UnconditionalWorkspace):
    '''compute conditional workspace.

    use only those parts of the workspace that contain both a segment
    and an annotation.
    '''

    is_conditional = True

    def __call__( self, segments, annotations, workspace ):
        
        # restrict workspace to those segments that contain both annotations and segments
        temp_workspace = workspace.clone()
        temp_workspace.filter( annotations )
        temp_workspace.filter( segments )

        return self.filter( segments, annotations, temp_workspace )

class ConditionalWorkspaceCentered(UnconditionalWorkspace):
    '''a workspace centered around segments/annotations.'''

    is_conditional = True

    def __init__( self, extension=None, expansion=None):
        self.extension = extension
        self.expansion = expansion
        if self.extension == self.expansion == None:
            raise ValueError( "need to specify either expansion or extension" )

    def __call__(self, segments, annotations, workspace):
        
        # build workspace from annotations
        temp_workspace = self.getCenter(segments, annotations).clone()
        if self.extension != None:
            temp_workspace.extend(self.extension)
        else:
            temp_workspace.expand(self.expansion)
        temp_workspace.normalize()

        # intersect with global workspace
        temp_workspace.intersect( workspace )

        return self.filter( segments, annotations, temp_workspace )

class ConditionalWorkspaceAnnotationCentered( ConditionalWorkspaceCentered ):
    '''compute conditional workspace.

    workspace is derived from the locations of the annotations.
    '''
    def getCenter( self, segments, annotations):
        return annotations

class ConditionalWorkspaceSegmentCentered( ConditionalWorkspaceCentered ):
    '''compute conditional workspace.

    workspace is derived from the locations of the segments
    '''
    # workspace needs to be computed on a per-segment basis
    is_conditional = False
    def getCenter( self, segments, annotations):
        return segments

############################################################
############################################################
############################################################
## compute counts
############################################################
def _append( counts, x,y,z): counts[y].append(z)
def _set( counts, x,y,z) : counts[x].__setitem__( y, z )
def _defdictfloat(): return collections.defaultdict(float)

def computeCounts(counter, 
                  aggregator, 
                  segments, 
                  annotations, 
                  workspace,
                  workspace_generator,
                  append = False):
    '''collect counts from *counter* between all combinations of
    *segments* and *annotations*.

    *aggregator* determines how values are combined across isochores.

    If *append* is set, values for each track in *segments* are
    appended to a list for each track in *annotations*. This is useful
    for samples.  Otherwise, a nested dictionary is returned.

    '''
        
    if append:
        counts = collections.defaultdict(list)
        f = _append
    else:
        counts = collections.defaultdict(_defdictfloat)
        f = _set

    isochores = workspace.keys()

    # collect counts per isochore and aggregate these
    for track in segments.tracks:
        segs = segments[track]
            
        for annotation in annotations.tracks:
            annos = annotations[annotation]

            temp_segs, temp_annos, temp_workspace = workspace_generator(
                segs, annos, workspace)
            vals = [counter(segs[isochore], annos[isochore], workspace[isochore])
                    for isochore in isochores ]
            f(counts, track, annotation, aggregator(vals))

    return counts

############################################################
############################################################
############################################################
## Quick parsing of tabular files (bed, gff)
############################################################
cdef class TupleProxy:
    '''Proxy class for access to parsed row as a tuple.

    This class represents a table row for fast read-access.
    '''

    cdef:
        char * data
        char ** fields
        int nfields
        int index
        int min_fields

    def __cinit__(self ):

        self.data = NULL
        self.fields = NULL
        self.index = 0
        self.nfields = 0
        self.min_fields = 0

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        if self.data != NULL:
            free(self.data)
            
        self.data = buffer
        self.update( buffer, nbytes )

    cdef present( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        if self.data != NULL:
            free(self.data)
            self.data = NULL

        self.update( buffer, nbytes )

    cdef copy(self, char * buffer, size_t nbytes):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        if self.data != NULL:
            free(self.data)

        cdef int s
        s = sizeof(char) * nbytes
        self.data = <char*>malloc( s )
        if not self.data:
            raise MemoryError("out of memory when allocation %i bytes" % sizeof(s))
        memcpy(<char*>self.data, buffer, s)
        self.update(self.data, nbytes)

    cdef update(self, char * buffer, size_t nbytes):
        '''update internal data.'''
        cdef char * pos
        cdef char * old_pos
        cdef int field
        cdef int max_fields
        field = 0

        if buffer[nbytes-1] != 0:
            raise ValueError( "incomplete line at %s: %s" % (buffer, buffer[nbytes-1]) )

        if self.fields != NULL:
            free(self.fields)

        max_fields = nbytes / 2
        self.fields = <char **>calloc( max_fields, sizeof(char *) )

        pos = buffer
        self.fields[0] = pos
        field += 1
        old_pos = pos

        while 1:

            pos = <char*>memchr( pos, '\t', nbytes )
            if pos == NULL: break
            pos[0] = '\0'
            pos += 1
            self.fields[field] = pos
            field += 1
            if field >= max_fields:
                raise ValueError("row too large - more than %i fields" % max_fields )
            nbytes -= pos - old_pos
            if nbytes < 0: break
            old_pos = pos

        self.nfields = field

        if self.nfields < self.min_fields:
            raise IOError("not enough fields, expected >%i, got %i" % \
                              (self.min_fields, self.nfields))

    def __getitem__(self, key):

        cdef int i
        i = key
        if i < 0: i += self.nfields
        if i >= self.nfields or i < 0:
            raise IndexError( "list index out of range" )
        return self.fields[i]

    def __len__(self):
        return self.nfields

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data)
        if self.fields != NULL:
            free(self.fields)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        """python version of next().
        """
        if self.index >= self.nfields:
            raise StopIteration
        self.index += 1
        return self.fields[self.index-1]

cdef class BedProxy(TupleProxy):

   cdef track
   def __init__(self, track ):
       self.track = track
       self.min_fields = 3

   property contig:
       def __get__(self):
           return force_str(self.fields[0])

   property start:
       def __get__(self):
           return atol(self.fields[1])

   property end:
       def __get__(self):
           return atol(self.fields[2])

   property name:
       def __get__(self):
           if self.nfields >= 4:
               return self.fields[3]
           else:
               return None

   property track:
       def __get__(self):
           return self.track

class Track(object):
    '''bed track information.'''
    def __init__(self, line):
        r= re.compile('([^\s=]+) *= *("[^"]*"|[^ ]*)')

        self._d = {}
        for k, v in r.findall(line[:-1]):
            if v[:1]=='"':
                self._d[k] = v[1:-1]
            else:
                self._d[k] = v

        self._line = line[:-1]

    def __str__(self):
        return self._line

    def __getitem__(self, key):
        return self._d[key]

    def __setitem__(self, key,val):
        self._d[key] = val


class tsv_iterator:
    '''iterate over ``infile``.

    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile):
        self.infile = infile

    def __iter__(self):
        return self

    def preparse(self, line):
        return True

    def create(self):
        return TupleProxy()

    def __next__(self):

        cdef char * b

        cdef TupleProxy r
        cdef size_t nbytes

        track = None

        while 1:

            line = self.infile.readline()
            if not line:
                break

            # skip comments
            if line.startswith("#") or line.startswith("\n"):
                continue

            if not self.preparse(line):
                continue

            # string conversion - b references internal string in python object

            # nbytes includes new line but not \0
            nbytes = len(line)

            line = force_bytes(line)
            b = line

            # make sure that entry is complete
            if b[nbytes-1] != '\n':
                raise ValueError( "incomplete line at %s" % line )
        
            # chop off newline
            b[nbytes-1] = '\0'

            # create a copy
            r = self.create()
            r.copy(b, nbytes)
            return r

        raise StopIteration

    # Python 2
    next = __next__
    

class bed_iterator(tsv_iterator):
    '''iterate over ``infile``.

    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile ):
        tsv_iterator.__init__(self, infile)
        self.track = None

    def preparse(self, line):
        if line.startswith("track"):
            self.track = Track(line)
            return False
        return True

    def create( self ):
        return BedProxy(self.track)


def readFromBed(filenames,
                bint allow_multiple=False,
                bint ignore_tracks=False):
    '''read SegmentLists from one or more bed files.

    Segment lists are grouped by `track` and `contig`..

    If there is no track information in the BED file, the name
    attribute (column 4) is taken instead. If the BED file has
    only three columns, the ``track`` will be called ``default``.

    Arguments
    ---------
    allow_multiple : bool
        If True, tracks can be spread over multiple files. The
        default is to raise an error if a track of the same name
        appears in more than one file.
    ignore_tracks : bool
        If True, ignore track information. There will only
        be a single track called "merged".

    '''
    cdef SegmentList l
    cdef BedProxy bed
    cdef Position lineno

    segment_lists = collections.defaultdict(IntervalDictionary)

    tracks = {}

    if type(filenames) == str:
        filenames = [filenames]

    for filename in filenames:
        with IOTools.openFile(filename, "r") as infile:
            default_name = os.path.basename(filename)
            lineno = 0
            try:
                for bed in bed_iterator(infile):
                    if ignore_tracks:
                        name = "merged"
                    else:
                        if bed.track:
                            try:
                                name = bed.track["name"]
                            except KeyError:
                                raise KeyError(
                                    "track without field 'name' in file '%s'" % filename)
                        elif bed.name: 
                            name = force_str(bed.name)
                        else:
                            name = default_name
                    
                    if name in tracks:
                        if tracks[name] != filename: 
                            if allow_multiple:
                                E.warn(
                                    "track '%s' in multiple filenames: %s and %s" %
                                    (name, tracks[name], filename))
                            else:
                                raise ValueError(
                                    "track '%s' in multiple filenames: %s and %s" %
                                    (name, tracks[name], filename))
                            tracks[name] = filename
                    else:
                        tracks[name] = filename

                    l = segment_lists[name][bed.contig]
                    l.add(atol(bed.fields[1]), atol(bed.fields[2]))
                    lineno += 1

            except IOError, msg:
                raise IOError("malformatted entry in line %s:%i, msg=%s" % (filename,
                                                                            lineno,
                                                                            msg))

    return segment_lists


cdef class IntervalContainer(object):
    '''generic container class representing a collection of SegmentList objects.
    
    The base class implements sending/retrieving data to/from shared memory.
    '''

    def __init__(self):
        self.name = None
        self.shared_fd = -1
        self.shared_fn = None

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def _getSharedMemoryName(self, key=None):
        '''return name of shared memory segment.

        If *key* is given that will be used as the suggested name
        for the memory segment. Otherwise, self.name will be used.
        
        Either way, the pid of the process will be added to the name.
        '''
        if key == None:
            assert self.name != None, \
                    "IntervalContainer: tried to build shared memory " \
                    "name without name defined"
            return "/" + str(os.getpid()) + "_" + self.name
        else:
            if key.startswith("/"): key = key[1:]
            return "/" + str(os.getpid()) + "_" + key

    def share(self, filename=None):
        '''setup collection for sharing.
        
        All contents are shared into a single memory mapped
        file.
        '''
        # determine size off memory required
        cdef off_t nbytes = self.counts() * sizeof(Segment)

        # open file in shared memory
        cdef int fd
        filename = self._getSharedMemoryName(filename)
        assert filename.startswith("/")
        f = force_bytes(filename)
        cdef char * c_filename = f
        fd = shm_open(c_filename,
                      O_CREAT | O_RDWR, 
                      S_IRUSR | S_IWUSR)
        if fd == -1:
            error = errno
            raise OSError("could not create shared memory at %s; "
                          "ERRNO=%i" % (filename, error ))

        # save file descriptor
        self.shared_fn = filename
        self.shared_fd = fd

        # resize file
        if ftruncate(fd, nbytes) == -1:
            error = errno
            raise OSError("could not resize memory at %s; ERRNO=%i" %
                          (filename, error))
        
        # setup memory map
        cdef void * mm
        mm = mmap(NULL, nbytes,
                  PROT_READ | PROT_WRITE, 
                  MAP_SHARED, 
                  fd, 
                  0);
        
        if mm == <void *>-1:
            raise ValueError("could not create memory mapped file")

        # copy all segment lists to shared memory
        cdef off_t offset = 0
        cdef CoordinateList clist
        for clist in self.getSegmentLists():
            offset = clist.toMMAP(mm, fd, offset)
        
    def __dealloc__(self):

        cdef int error
        cdef CoordinateList clist

        if self.shared_fd != -1:
            # if IntervalCollection is destroyed that has
            # used memory mapped SegmentLists, the 
            # SegmentLists are moved back to private
            # memory.
            # This is potentially unnecessary copying
            # if all SegmentLists are subsequently destroyed
            # but necessary if there the SegmentLists are
            # referenced outside the IntervalCollection.
            for clist in self.getSegmentLists():
                clist.fromMMAP()

            # disconnect from mmap
            # munmap( self.mmap, 
            #        self.mmap_bytes)
            fn = force_bytes(self.shared_fn)
            fd = shm_unlink(fn)
            error = errno
            if fd == -1:
                raise OSError(
                    "IntervalCollection.__dealloc__(): could not unlink shared memory "
                    "%s: ERRNO=%i" % (self.shared_fn, error))
            else:
                E.debug("IntervalCollection.__dealloc__(): freeing shared memory "
                        "at %s" % self.shared_fn)
                
    def unshare(self):
        '''revert sharing.'''
        cdef CoordinateList clist

        assert self.shared_fd != -1, "unsharing shared IntervalCollection"

        for clist in self.getSegmentLists():
            clist.fromMMAP()
                
        # disconnect from mmap
        # munmap( self.mmap, 
        #        self.mmap_bytes)
                
        fn = self.shared_fn 
        fd = shm_unlink(fn)
        error = errno
        if fd == -1:
            raise OSError(
                "IntervalCollection.unshare(): could not unlink shared memory "
                "%s: ERRNO=%i" % (self.shared_fn,error))
        else:
            E.debug("IntervalCollection.unshare(): freeing shared memory "
                    "at %s" % self.shared_fn)
        
        self.shared_fd = -1

    def sum(self):
        '''return sum of all segment lists.'''
        s = 0
        for segmentlist in self.getSegmentLists():
            s += segmentlist.sum()
        return s

    def counts(self):
        '''return number of all segments in all segments lists.'''
        s = 0
        for segmentlist in self.getSegmentLists():
            s += len(segmentlist)
        return s

    def sort(self):
        '''sort all intervals lists.'''
        for segmentlist in self.getSegmentLists():
            segmentlist.sort()

    def extend(self, extension):
        '''extend each interval by a certain amount.'''
        cdef SegmentList segmentlist
        for segmentlist in self.getSegmentLists():
            segmentlist.extend_segments(extension)

    def expand(self, expansion):
        '''expand each interval by a certain amount.'''
        cdef SegmentList segmentlist
        for segmentlist in self.getSegmentLists():
            segmentlist.expand_segments(expansion)

    def normalize(self):
        '''normalize all segment lists.'''
        for segmentlist in self.getSegmentLists():
            segmentlist.normalize()

    def check(self):
        '''check all intervals lists.'''
        self.sort()


cdef class IntervalDictionary(IntervalContainer):
    '''a collection of intervals.
    
    Intervals (objects of type :class:`SegmentList`) are organized
    in a hierarchical dictionary by isochore.
    '''

    def __init__(self, name = None, unreduce = None ):
        self.intervals = collections.defaultdict(SegmentList)
        self.shared_fd = -1
        self.name = name
        self.shared_fn = None
        if unreduce:
            (self.name, self.intervals) = unreduce
            # shared_fd and shared_fn remain unset - marks object as slave
            
    def __reduce__(self):
        return (buildIntervalDictionary, (self.name, 
                                          self.intervals ))


    def getSegmentLists( self ):
        '''yield all segment lists.'''
        for contig, segmentlist in self.intervals.items():
            yield segmentlist

    def intersect(self, other):
        '''intersect with intervals in other.'''
        for contig, segmentlist in self.intervals.items():
            if contig in other:
                segmentlist.intersect(other[contig])
            else:
                del self.intervals[contig]

    def truncate( self, other ):
        '''truncate intervals with intervals in other.'''
        for contig, segmentlist in self.intervals.items():
            if contig in other:
                segmentlist.truncate(other[contig])
            else:
                del self.intervals[contig]

    def __len__(self):
        return len(self.intervals)

    def __str__(self):
        return ";".join(["%s:%i,%i" % (x, len(y), y.sum())
                         for x,y in self.intervals.items()])

    def __delitem__(self,key):
        del self.intervals[key]

    def keys(self): return self.intervals.keys()

    def add( self, contig, segmentlist):
        self.intervals[contig] = segmentlist

    def items(self ):
        return self.intervals.items()

    def __getitem__(self, key):
        return self.intervals[key]

    def __setitem__(self, key, val):
        self.intervals[key] = val

    def __contains__(self, key ):
        return key in self.intervals

    def items(self):
        return self.intervals.items()

    def clone(self):
        '''return a copy of the data.'''
        r = IntervalDictionary()
        for contig, segmentlist in self.intervals.items():
            r[contig] = segmentlist.clone()
        return r

    def filter(self, other):
        '''remove all intervals not overlapping with intervals in other.'''
        for contig, segmentlist in self.intervals.items():
            if contig in other:
                segmentlist.filter(other[contig])
            else:
                del self.intervals[contig]

    def toIsochores(self, isochores, truncate=False):
        '''split per-contig segmentlists into per-isochore segmentlist.

        If *truncate* is given, the segments are truncated at isochore
        boundaries.

        The IntervalDictionary is modified in-place.
        '''
        for contig in list(self.intervals.keys()):
            segmentlist = self.intervals[contig]
            for other_track, other_vv in isochores.items():
                newlist = segmentlist.clone()
                if truncate:
                    newlist.intersect(other_vv[contig])
                else:
                    newlist.filter(other_vv[contig])
                isochore = "%s.%s" % (contig, other_track)
                self.intervals[isochore] = newlist
            del self.intervals[contig]

    def fromIsochores(self):
        '''merge isochores into contigs'''
        new = collections.defaultdict(SegmentList)
        normalize = False
        # isochores might or might not be present
        for isochore, segmentlist in self.intervals.items():
            isochore = isochore.strip()
            if "." in isochore and isochore != ".":
                contig, iso = isochore.split(".")
                new[contig].extend(segmentlist)
                normalize = True
            else:
                new[isochore] = segmentlist 

        if normalize:
            # merge adjacent intervals (and normalize)
            for x in new.values(): 
                x.merge(0)
                
        self.intervals = new

def buildIntervalDictionary(*args):
    '''unpickling - return a rebuilt SamplerShift object.'''
    return IntervalCollection(unreduce=args)

#####################################################################
#####################################################################
#####################################################################
## 
#####################################################################
cdef class IntervalCollection(IntervalContainer):
    '''a collection of intervals.

    Intervals (objects of type :class:`SegmentList`) are organized
    in hierarchical dictionary first by track and then by isochore.
    '''


    def __init__(self, name=None, unreduce=None):
        self.intervals = collections.defaultdict(IntervalDictionary)
        self.name = name
        self.shared_fd = -1
        self.shared_fn = None
        if unreduce:
            (self.name, self.intervals) = unreduce
            # shared_fd and shared_fn remain unset - marks object as slave

    def __reduce__(self):
        return (buildIntervalCollection, (self.name, 
                                          self.intervals ))

    def setName(self, name):
        """set name of collection."""
        self.name = name

    def getSegmentLists(self):
        '''yield all segment lists.'''
        for track,v in self.intervals.items():
            for contig, segmentlist in v.items():
                yield segmentlist

    def load(self, filenames, allow_multiple=False, ignore_tracks=False):
        '''load segments from filenames.'''
        self.intervals = readFromBed(filenames,
                                     allow_multiple=allow_multiple,
                                     ignore_tracks=ignore_tracks)

    def save(self, outfile, prefix = "", **kwargs):
        '''save in bed format to *outfile*.
        
        Each interval set will be saved as a different track with optional
        prefix *prefix*.

        Additional *kwargs* will be added to the tracks line as key=value pairs.
        '''

        for track, vv in self.intervals.items():
            outfile.write(
                "track name=%s%s %s\n" %
                (prefix, track,
                 " ".join(["%s=%s" % (x,y) for x,y in kwargs.items()])))
            for contig, segmentlist in vv.items():
                for start, end in segmentlist:
                    outfile.write("%s\t%i\t%i\n" % (contig, start, end))

    def normalize(self):
        '''normalize segment lists individually.

        Remove empty contigs.
        '''

        for track, vv in self.intervals.items():
            for contig in vv.keys():
                segmentlist = vv[contig]
                if len(segmentlist) == 0: 
                    del vv[contig]
                else:
                    segmentlist.normalize()

    def outputStats(self, outfile):
        '''output segment statistics.'''

        outfile.write( "section\ttrack\tcontig\tnsegments\tlength\n" )

        cdef Position total_length = 0
        cdef Position total_segments = 0
        cdef Position length, segments
        for track, vv in self.intervals.items():
            total_length, total_segments = 0, 0
            for contig, segmentlist in vv.items():
                segments = len(segmentlist)
                length = segmentlist.sum()
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % segments,
                         "%i" % length ) ) + "\n" )
                total_length += length
                total_segments += segments
            outfile.write("\t".join( \
                (self.name, track, "total", 
                         "%i" % total_segments,
                         "%i" % total_length ) ) + "\n" )
                 
                
    def merge(self, delete=False):
        '''merge all tracks into a single segment list
        creating a new track called 'merged'.

        Overlapping intervals will be merged.
        
        .. note::

            The merged track will not be sorted or
            normalized.

        Arguments
        ----------
        delete : bool
            If True, all tracks but the merged track
            will be deleted.

        '''
        merged = IntervalDictionary()
        
        for track in self.intervals.keys():
            vv = self.intervals[track]
            for contig, segmentlist in vv.items():
                merged[contig].extend(segmentlist)
            if delete:
                del self.intervals[track]

        self.intervals["merged"] = merged

    def collapse( self ):
        '''collapse all tracks into a single segment list
        creating a new track called 'collapsed'.

        The collapsed track contains the intersection of
        all segment lists.
        '''
        result = IntervalDictionary()
        
        # get list of contigs present in all data sets
        contigs = collections.defaultdict( int )
        for track, vv in self.intervals.items():
            for contig in vv.keys():
                contigs[contig] += 1

        ntracks = len( self.intervals )
        shared_contigs = set( [ x for x,y in contigs.items() if y == ntracks ] )
        
        for track, vv in self.intervals.items():
            for contig, segmentlist in vv.items():
                if contig not in shared_contigs: continue
                if contig not in result:
                    result[contig] = segmentlist.clone()
                else:
                    result[contig].intersect(segmentlist)
                    
        self.intervals["collapsed"] = result

    def countsPerTrack(self):
        '''return number of all segments in all segments lists.'''
        counts = {}
        for track, vv in self.intervals.items():
            s = 0
            for contig, segmentlist in vv.items():
                s += len(segmentlist)
            counts[track] = s 
        return counts
        
    def intersect(self, other):
        '''intersect with intervals in other.'''
        for track, vv in self.intervals.items():
            for contig, segmentlist in vv.items():
                if contig in other:
                    segmentlist.intersect( other[contig] )
                else:
                    del vv[contig]

    def filter(self, other):
        '''remove all intervals not overlapping with intervals in other.'''
        for track, vv in self.intervals.items():
            for contig, segmentlist in vv.items():
                if contig in other:
                    segmentlist.filter( other[contig] )
                else:
                    del vv[contig]

    def restrict(self, restrict):
        '''remove all tracks except those in restrict.'''
        if restrict in (list, tuple, set):
            r = set(restrict)
        else:
            r = set([restrict,])

        keys = list(self.intervals.keys())
        for track in keys:
            if track not in r:
                del self.intervals[track]

    def toIsochores(self, isochores, truncate=False):
        '''split per-contig segmentlists into per-isochore segmentlist.

        If *truncate* is given, the segments are truncated at isochore
        boundaries.

        The IntervalCollection is modified in-place.
        '''
        for track, vv in self.intervals.items():
            vv.toIsochores(isochores, truncate)

    def fromIsochores(self):
        '''merge isochores together.'''
        for track, vv in self.intervals.items():
            vv.fromIsochores()

    def toPositions(self, method="mid-point"):
        """convert SegmentLists to PositionLists"""
        for track, vv in self.intervals.items():
            for contig, segmentlist in vv.items():
                p = PositionList()
                p.fromSegmentList(segmentlist, method=method)
                vv[contig] = p

    def clone(self):
        '''return a copy of self.'''
        new = IntervalCollection(self.name)
        for track,v in self.intervals.items():
            for contig, segmentlist in v.items():
                new.add(track, contig, segmentlist.clone())
        return new

    @property
    def tracks(self):
        return self.intervals.keys()

    def __len__(self):
        return len(self.intervals)

    def __delitem__(self,key):
        del self.intervals[key]

    def keys(self):
        return self.intervals.keys()

    def add(self, track, contig, segmentlist):
        self.intervals[track][contig] = segmentlist

    def __getitem__(self, key):
        return self.intervals[key]

    def __contains__(self, key):
        return key in self.intervals

    def __str__(self):
        return "%s:%s" % (
            self.name,
            ",".join( ["%s:%s" % (x,y) for x,y in self.intervals.items()]))

    def items(self):
        return self.intervals.items()
    
    def items(self):
        return self.intervals.items()
    
    def outputOverlapStats(self, outfile, other):

        outfile.write("section\ttrack\tcontig\toverlap\tlength\tdensity\n")
        
        for track, vv in self.intervals.items():
            for contig, segmentlist in vv.items():
                length = segmentlist.sum()
                if length == 0: continue
                overlap = segmentlist.overlapWithSegments(other[contig])
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % overlap,
                         "%i" % length,
                         "%f" % (float(overlap) / length)) ) + "\n" )

def buildIntervalCollection(*args):
    '''unpickling - return a rebuilt SamplerShift object.'''
    return IntervalCollection(unreduce=args)


cdef class Samples(object):
    '''a collection of samples.

    Samples :class:`IntervalCollections` identified by track and sample_id.
    '''
    cdef dict samples

    def __init__(self ):
        '''create a new SampleCollection.

        If cache is given, samples will be stored persistently on disk.
        '''
        self.samples = {}

    def add( self, track, sample_id, isochore, segmentlist ):
        '''add a new *sample* for *track* and *isochore*, giving it *sample_id*.'''
        if track not in self.samples:
            self.samples[track] = IntervalCollection( track )
        self.samples[track].add( sample_id, isochore, segmentlist )

    def hasSample( self, track, sample_id, isochore ):
        '''return true if cache has sample.'''
        if track not in self.samples: return False
        if sample_id not in self.samples[track]: return False
        return isochore in self.samples[track][sample_id]

    def __contains__(self, track ):
        return track in self.samples

    def __getitem__(self, track ):
        '''return all samples for track (as an :class:`IntervalCollection`)'''
        return self.samples[track]

    def load( self, track, sample_id, isochore ):
        raise ValueError( "loading called for uncached data" )

    def __delitem__(self, key ):
        E.debug("releasing memory for sample %s" % key )
        del self.samples[key]

    def __len__(self):
        return len(self.samples)

cdef class SamplesFile( Samples ):
    '''a collection of samples.

    Samples :class:`IntervalCollections` identified by track and sample_id.
    '''
    def __init__(self, filenames, regex ):
        '''create a new SampleCollection.
        
        Samples are read from *filenames* at startup.
        The track name is given by applying the regular
        expression to the pattern.
        '''
        Samples.__init__(self)
        
        for filename in filenames:
            track = regex.search(filename).groups()[0]
            s = IntervalCollection( track )
            s.load( filename)
            self.samples[track] = s
            
    def load(self, track, sample_id, isochore ):
        return True


cdef class SamplesCached( Samples ):
    '''a collection of samples.

    Samples :class:`IntervalCollections` identified by track and sample_id.
    '''
    cdef FILE * fcache
    cdef FILE * findex
    cdef dict index
    cdef char * filename

    def __init__(self, filename ):
        '''create a new SampleCollection.

        If cache is given, samples will be stored persistently on disk.
        '''
        Samples.__init__(self)
        
        self.filename = filename
        tmp = self.filename + ".idx"

        if not os.path.exists( filename ):
            self.fcache = fopen( filename, "wb" )
            self.findex = fopen( tmp, "wb" )
            self.index = {}
        else:
            self.fcache = fopen( filename, "rb" )
            self.loadIndex()
            self.findex = fopen( tmp, "rb" )

    def loadIndex( self ):
        '''load index from cache.

        The index is loaded from a new file.
        '''
        E.debug( "loading index from %s" % self.filename )
        cdef char * ckey
        cdef off_t pos
        cdef char keylen
        cdef FILE * findex
        self.index = {}
        tmp = self.filename + ".idx"
        findex = fopen( tmp, "rb" )

        ckey = <char*>calloc(sizeof(char) * 256, 1)
        x = 0
        while not feof(findex):
            fread(&keylen, sizeof( char), 1, findex)
            if feof(findex):
                break
            fread(ckey, sizeof( char), keylen, findex)
            fread(&pos, sizeof(off_t), 1, findex)
            # converts to persistent python object
            self.index[ckey] = pos
            
        fclose(findex)
        free(ckey)

        E.debug( "loaded index from %s: %i items" % (self.filename,len(self.index)))

    def add( self, track, sample_id, isochore, segmentlist):
        '''add a new *sample* for *track* and *isochore*, giving it *sample_id*.'''

        Samples.add( self, track, sample_id, isochore, segmentlist )
        cdef SegmentList seglist
        seglist = segmentlist
        cdef off_t pos
        cdef char * key
        cdef char keylen
        
        l = len(seglist)
        if l == 0: return

        pos = ftell(self.fcache)

        # cache structure is:
        # 1 * sizeof(unsigned char) - key length (max 255 chars)
        # 1 * sizeof(Position) - number of segments (nsegments)
        # nsegments * sizeof(Segment) - the segment list
        tempkey = self.toKey(track, sample_id, isochore)
        assert len(tempkey) <= 255
        key = tempkey
        keylen = strlen(key) + 1

        self.index[key] = pos

        fwrite(&seglist.nsegments, sizeof(Position), 1, self.fcache)
        toCompressedFile(<unsigned char *>seglist.segments,
                         sizeof(Segment) * seglist.nsegments,
                         self.fcache)

        # write to index
        fwrite(&keylen, sizeof(char), 1, self.findex)
        fwrite(key, sizeof(char), keylen, self.findex)
        fwrite(&pos, sizeof(off_t), 1, self.findex)

    def toKey( self, track, sample_id, isochore ):
        return "%s-%s-%s" % (track,sample_id,isochore) 

    def hasSample(self, track, sample_id, isochore):
        '''return true if cache has sample.'''
        return self.toKey(track, sample_id, isochore) in self.index

    def __dealloc__(self):
        fclose( self.fcache)
        fclose( self.findex)

    def load(self, track, sample_id, isochore):
        '''load data into memory'''
        cdef off_t pos
        cdef SegmentList seglist
        cdef Position nsegments

        tempkey = self.toKey(track, sample_id, isochore)
        pos = self.index[tempkey]
        fseek(self.fcache, pos, SEEK_SET)
        fread(&nsegments, sizeof(Position), 1, self.fcache)
        seglist = SegmentList(allocate=nsegments)
        seglist.nsegments = nsegments
        fromCompressedFile(<unsigned char*> seglist.segments,
                           sizeof(Segment) * seglist.nsegments,
                           self.fcache)

        Samples.add(self, track, sample_id, isochore, seglist)

############################################################
############################################################
############################################################
## FDR computation - obsolete
############################################################
cdef double computeFalsePositiveRate( EnrichmentStatistics ** allstats,
                                      Position nresults,
                                      Position nsamples,
                                      double pvalue ):
    '''return the number of expected number of false positives
    for less than or equal to *pvalue*.
    '''

    cdef Position y, nsample, nfp, total_nfp
    cdef double efp

    total_nfp = 0

    # Compute for each sample the number of simulated
    # data points that are called significant at pvalue.
    # As there are the same number of "terms" in each sample
    # that are tested, a simple summation suffices.
    for nsample from 0 <= nsample < nsamples:
        nfp = 0
        for y from 0 <= y < nresults:
            if isSampleSignificantAtPvalue( allstats[y], nsample, pvalue ): 
                nfp += 1
        total_nfp += nfp

    efp = <double>total_nfp / nsamples

    return efp

cpdef computeFDR(annotator_results):
    '''compute an experimental fdr across all segments and annotations.

    The experimental fdr is given by

    E(FP) = average number of terms in each sample run with P-Value <= p
        aka: experimental number of false positives

    R = number of nodes in observed data, that have a P-Value of less than or 
        equal to p.
        aka: pos=positives in observed data

    fdr = E(FP)/R

    The results are added to annotator_results.
    '''
    fdr_cache = {}

    cdef Position nresults
    cdef AnnotatorResult r
    cdef EnrichmentStatistics * r1
    cdef EnrichmentStatistics ** allstats
    cdef Position nsample, nfp, R, x,y, total_nfp
    cdef Position * nfps

    cdef double pvalue, efp

    nresults = len(annotator_results)

    # collect results
    allstats = <EnrichmentStatistics **>calloc( nresults, sizeof(EnrichmentStatistics *))

    for x from 0 <= x < nresults:
        r = annotator_results[x]
        allstats[x] = r.stats

    all_pvalues = numpy.array( [ r.stats.pvalue for r in annotator_results ], dtype = numpy.float )
    all_pvalues.sort()

    for x from 0 <= x < nresults:

        E.debug( "progress: %i/%i " % (x+1, nresults))
        r1 = allstats[x]

        pvalue = r1.pvalue
        
        if pvalue in fdr_cache:
            r1.qvalue = fdr_cache[ pvalue ]
            continue

        efp = computeFalsePositiveRate( allstats, 
                                        nresults,
                                        r1.nsamples,
                                        pvalue )

        # number of positives at P-Value
        R = numpy.searchsorted( all_pvalues, r1.pvalue )
        while R < nresults and all_pvalues[R] <= pvalue:
            R += 1

        r1.qvalue = dmin( 1.0, dmax( 1.0 / r1.nsamples, efp / R))
        # print "pvalue=", r1.pvalue, "qvalue=", r1.qvalue, "efp=", efp, "nfp=", total_nfp, "nsamples=", \
            # r1.nsamples, "R=", R, "vals=", all_pvalues[max(0,R-3):R+3]

        fdr_cache[pvalue] = r1.qvalue
        
        break

@cython.profile(False)
cdef inline int isSampleSignificantAtPvalue( EnrichmentStatistics * stats, 
                                             Position sample_id, 
                                             double pvalue ):
    '''return True, if sample sample_id would be called
    significant at threshold *pvalue*

    This method is fast for small pvalues, but slow for large
    pvalues because the method does not a full search.

    This method works by scanning the first/last samples
    until one is found that is larger/smaller than the
    value of samples[sample_id].
    '''
    cdef double pval, min_pval, val
    cdef Position l

    l = stats.nsamples
    min_pval = 1.0 / l

    cdef int idx
    idx = stats.sample2sorted[sample_id]
    val = stats.samples[sample_id]
    if val > stats.expected:
        # over-representation
        while idx > 0 and stats.samples[stats.sorted2sample[idx]] == val: 
            idx -= 1
        idx = l - (idx + 1)
    elif val < stats.expected:
        # under-representation
        while idx < l and stats.samples[stats.sorted2sample[idx]] == val: 
            idx += 1

    pval = float(idx) / l

    # = is important, such that a sample itself is called significant
    # print "val=", val, "idx=", idx, "sample_id=", sample_id, "pval=", pval, "refpval=", pvalue, "true=", pval <= pvalue
    return dmax( min_pval, pval ) <= pvalue


# import tables
# import warnings
# warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)


# class SamplesPytables( object ):
#     '''a collection of samples.

#     Samples :class:`IntervalCollections` identified by track and sample_id.
#     '''

#     def __init__(self, cache = None ):
#         '''create a new SampleCollection.

#         If cache is given, samples will be stored persistently on disk.
#         '''
#         if cache:
#             self.cache = tables.openFile( cache, mode = "a", title = "Sample cache")
#             self.filters = tables.Filters(complevel=5, complib='zlib')
#         else:
#             self.cache = None

#         self.samples = collections.defaultdict( IntervalCollection )

#     def add( self, track, sample_id, isochore, sample ):
#         '''add a new *sample* for *track* and *isochore*, giving it *sample_id*.'''
#         if track not in self.samples:
#             self.samples[track] = IntervalCollection( track )

#         self.samples[track].add( sample_id, isochore, sample )

#         if self.cache:
#             l = len(sample)
#             if l == 0: return

#             try:
#                 loc = self.cache.getNode( "/%s/%i" % (track, sample_id) )
#             except tables.exceptions.NoSuchNodeError:
#                 loc = self.cache.createGroup( "/%s/%i" % (track, sample_id),
#                                               "%s-%i" % (track, sample_id),
#                                               "%s-%i" % (track, sample_id),
#                                               createparents = True )

#             carr = self.cache.createCArray( loc,
#                                             isochore,
#                                             tables.UInt32Atom(),
#                                             shape=( l, 2),
#                                             filters = self.filters )

#             for x, c in enumerate( sample ):
#                 carr[x] = [c[0], c[1]]

#             carr.flush()

#     def hasSample( self, track, sample_id, isochore ):
#         '''return true if cache has sample.'''
#         if self.cache:
#             return "/%s/%i/%s" % (track, sample_id,isochore) in self.cache
#         else:
#             if track not in self.samples: return False
#             if sample_id not in self.samples[track]: return False
#             return isochore in self.samples[track][sample_id]

#     def save( self ):
#         '''save full interval collection in cache.'''

#         if not self.cache: return

#         for track, samples in self.samples.items():
#             group = self.cache.createGroup( self.cache.root, track, track)
#             for sample_id, sample in samples.items():
#                 subgroup = self.cache.createGroup( group, str(sample_id), str(sample_id) )
#                 for isochore, seglist in sample.items():
#                     l = len(seglist)
#                     if l == 0: continue
#                     carr = self.cache.createCArray( subgroup,
#                                                 isochore,
#                                                 tables.UInt32Atom(),
#                                                 shape=( l, 2),
#                                                 filters = self.filters )

#                     for x, c in enumerate( seglist ):
#                         carr[x] = [c[0], c[1]]
#                     carr.flush()




#     def __del__(self):
#         if self.cache:
#             self.cache.close()

#     def __getitem__(self, track ):
#         '''return all samples for track (as an :class:`IntervalCollection`)'''
#         return self.samples[track]

