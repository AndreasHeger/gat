# cython: embedsignature=True
# cython: profile=False

import types, collections, re, os, math, random

cimport cython

#####################################################
#####################################################
#####################################################
## numpy import
## both import and cimport are necessary
#####################################################
import numpy
cimport numpy
DTYPE_INT = numpy.int
ctypedef numpy.int_t DTYPE_INT_t
DTYPE_FLOAT = numpy.float
ctypedef numpy.float_t DTYPE_FLOAT_t

def getSegmentSize():
    '''return size of coordinate and size of a segment.'''
    return sizeof(Segment) // 2, sizeof(Segment)

# min/max are not optimized, so declare them as C functions
# declare as signed comparisons as a Position might be negative
@cython.profile(False)
cdef inline PositionDifference lmin( PositionDifference a, PositionDifference b):
    if a < b: return a
    return b
@cython.profile(False)
cdef inline PositionDifference lmax( PositionDifference a, PositionDifference b):
    if a > b: return a
    return b

@cython.profile(False)
cdef inline double dmin( double a, double b):
    if a < b: return a
    return b

@cython.profile(False)
cdef inline double dmax( double a, double b):
    if a > b: return a
    return b

@cython.profile(False)
cdef inline PositionDifference segment_overlap( Segment a, Segment b ):
    return lmax(0, <PositionDifference>lmin( a.end, b.end) - <PositionDifference>lmax(a.start, b.start))

@cython.profile(False)
cdef inline PositionDifference range_overlap( Position astart, Position aend, Position bstart, Position bend ):
    return lmax(0, 
                <PositionDifference>lmin( aend, bend) - 
                <PositionDifference>lmax(astart, bstart))

@cython.profile(False)
cdef inline PositionDifference segment_overlap_raw(  Segment a, Segment b ):
    return <PositionDifference>lmin(a.end, b.end) - <PositionDifference>lmax(a.start, b.start)

@cython.profile(False)
cdef inline PositionDifference segment_length(  Segment a ):
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
    return <int>((<double*>s1)[0] - (<double*>s2)[0])

cdef class SegmentList:
    '''list of segments.

    A list of (non-overlapping) segments.

    The intersection algorithms assume that segments
    are non-overlapping and sorted by start.

    Call normalize to merge overlapping segments.
    '''
    def __init__(self, 
                 int allocate = 0,
                 SegmentList clone = None,
                 iter = None,
                 normalize = False ):
        '''create empty list of segments.

        
        Creating is a lazy operation, no memory is allocated unless
        explitely given:

        *allocate* - empty list, allocate slots
        *clone* - copy an existing list
        *iter* - initialize list from iterator

        If *normalize* is set, the list will be normalized.
        '''
        cdef long idx, nsegments

        self.segments = NULL
        self.nsegments = 0
        self.allocated = 0
        # an empty list is normalized
        self.is_normalized = 1
        self.chunk_size = 10000

        if clone != None:
            self.nsegments = self.allocated = clone.nsegments
            self.segments = <Segment*>calloc( clone.nsegments, sizeof( Segment ) )
            memcpy( self.segments,
                    clone.segments,
                    clone.nsegments * sizeof(Segment ) )
            self.is_normalized = clone.is_normalized
        elif iter:
            a = tuple(iter)
            nsegments = len(a)
            self.nsegments = self.allocated = nsegments
            self.segments = <Segment*>calloc( self.nsegments, sizeof( Segment ) )
            idx = 0
            for start, end in a:
                self.segments[idx] = Segment( start, end )
                idx += 1
            self.is_normalized = 0
        elif allocate:
            self.allocated = allocate
            self.segments = <Segment*>calloc( allocate, sizeof( Segment ) )

        if normalize: self.normalize()

    cpdef sort( self ):
        '''sort segments.'''
        if self.nsegments == 0: return

        qsort( <void*>self.segments,
               self.nsegments,
               sizeof( Segment ),
               &cmpSegments )

    cpdef SegmentList extend( self, SegmentList other ):
        '''extend a list of segments with segments in other list.

        The list will not be normalized automatically - call
        :meth:`SegmentList.normalize`.
        '''
        cdef size_t new_size
        new_size = other.nsegments + self.nsegments
        if self.allocated == 0:
            self.segments = <Segment*>malloc( new_size * sizeof( Segment ) )
            assert self.segments != NULL
            self.allocated = new_size
        elif new_size >= self.allocated:
            self.segments = <Segment*>realloc( <void *>self.segments, new_size * sizeof( Segment ) )
            assert self.segments != NULL
            self.allocated = new_size
        memcpy( &self.segments[self.nsegments],
                 other.segments,
                 other.nsegments * sizeof(Segment ) )
        self.nsegments = new_size
        self.is_normalized = 0
        return self

    cdef _add( self, Segment segment ):
        '''add a new segment.

        The list will not be normalized automatically - call
        :meth:`SegmentList.normalize`.
        '''

        if self.allocated == 0:
            self.segments = <Segment*>malloc( self.chunk_size * sizeof( Segment ) )
            self.allocated = self.chunk_size
        elif self.nsegments == self.allocated:
            self.allocated *= 2
            self.segments = <Segment*>realloc( self.segments,
                                               self.allocated * sizeof( Segment ) )
        assert self.segments != NULL, "memory error"

        self.segments[self.nsegments] = segment
        self.nsegments += 1
        self.is_normalized = 0

    cpdef add( self, Position start, Position end ):
        cdef Segment segment
        assert start <= end, "attempting to add invalid segment %i-%i" % (start, end)
        segment = Segment( start, end)
        self._add( segment )

    cpdef trim_ends( self, Position pos, Position size, int forward ):
        '''trim segment list by removing *size* nucleotides from
        the segment that includes *pos*.

        The flag *forward* gives the direction.
        '''
        if self.nsegments == 0: return

        assert self.is_normalized, "trimming in non-normalized list"

        cdef int idx
        cdef Position l
        cdef PositionDifference s = size
        cdef Segment other, seg

        assert self.sum() > s, "trimming more than the total length (%i < %i)" % (self.sum(), s)
        #print "at start=%i, removing %i" % (self.sum(), s)

        other = Segment( pos, pos + 1)
        idx = self._getInsertionPoint( other )

        # wrap around
        if idx == self.nsegments: idx = 0
        if idx < 0: idx = self.nsegments - 1

        if forward:
            while s > 0:
                seg = self.segments[idx]
                l = segment_length( seg )
                if segment_length( seg ) < s:
                    self.segments[idx] = Segment(0,0)
                    s -= l
                else:
                    self.segments[idx] = Segment( seg.start + s, seg.end )
                    s = 0

                idx += 1
                if idx == self.nsegments: idx = 0
        else:
            while s > 0:
                seg = self.segments[idx]
                l = segment_length( seg )
                if segment_length( seg ) < s:
                    self.segments[idx] = Segment(0,0)
                    s -= l
                else:
                    self.segments[idx] = Segment( seg.start, <PositionDifference>seg.end - s )
                    s = 0

                idx -= 1
                if idx < 0: idx = self.nsegments - 1

        #print "at end=%i, removing %i" % (self.sum(), s)

    cdef _resize( self, int nsegments ):
        '''resize segment list to nsegments.'''
        # do not free all if there are no segments
        if nsegments == 0: nsegments = 1

        assert nsegments >= self.nsegments, "resizing will loose segments"
        
        if self.allocated == 0:
            self.segments = <Segment*>malloc( nsegments * sizeof( Segment ) )
        else:
            self.segments = <Segment*>realloc( self.segments,
                                               nsegments * sizeof( Segment ) )
        assert self.segments != NULL
        self.allocated = nsegments

    cdef insert( self, int idx, Segment seg ):
        '''insert Segment *seg* at position *idx*'''
        if idx < 0: raise ValueError( "only positive indices accepted (%i)" % idx )
        
        idx = max( idx, self.nsegments )
        if self.allocated == self.nsegments:
            self._resize( self.allocated * 2)

        memmove( &self.segments[idx+1],
                 &self.segments[idx],
                 sizeof( Segment ) * (self.nsegments - idx) )

        self.nsegments += 1
        self.segments[idx] = seg

    cpdef trim( self, Position pos, Position size ):
        '''trim segment list by removing *size* nucleotides
        starting at pos.

        If (pos, pos+size) is fully within a segment, this segment
        will be split.

        '''
        if self.nsegments == 0: return

        assert self.is_normalized, "trimming in non-normalized list"

        cdef int idx
        cdef PositionDifference l
        cdef PositionDifference s = size
        cdef PositionDifference p = pos
        cdef Segment other, seg

        assert self.sum() > s, "trimming more than the total length (%i < %i)" % (self.sum(), s)

        p = lmax( p, self.segments[0].start )
        p = lmin( p, self.segments[self.nsegments-1].end )

        other = Segment( p, p + 1)
        idx = self._getInsertionPoint( other )
        # wrap around
        if idx == self.nsegments:
            idx = 0
            p = self.segments[0].start
        
        #assert p >= self.segments[idx].start, "pos=%i, %i, %s, %s" % (p, idx, str(self.segments[idx]), str(self.asList()))
        #assert p < self.segments[idx].end, "pos=%i, %i, %s, %s" % (p, idx, str(self.segments[idx]), str(self.asList()))

        while s > 0:
            seg = self.segments[idx]
            
            l = <PositionDifference>self.segments[idx].end - p

            if s < l: 
                if p == self.segments[idx].start:
                    # truncate from left and stop
                    self.segments[idx].start = p + s
                    break
                else:
                    # remove from middle and stop
                    self.insert( idx + 1, Segment( p + s, self.segments[idx].end) )
                    self.segments[idx].end = p
                    break
            else:
                self.segments[idx].end = p
                s -= l

            idx += 1

            # wrap around
            if idx == self.nsegments: idx = 0
            p = self.segments[idx].start

        #print "at end=%i, removing %i" % (self.sum(), s)

    cpdef normalize( self ):
        '''merge all overlapping segments and remove empty segments.

        Adjacent segments are not merged.

        This function works in-place.
        '''

        cdef int idx, insertion_idx
        cdef Position max_end
        if self.nsegments == 0:
            self.is_normalized = 1
            return

        self.sort()

        insertion_idx = 0

        # skip to first non-empty segment
        idx = 0
        while idx < self.nsegments and self.segments[idx].start == self.segments[idx].end: 
            idx += 1

        # only empty segments - empty list, but declare normalized
        if idx == self.nsegments:
            self.nsegments = 0
            self.is_normalized = 1
            return

        self.segments[insertion_idx].start = self.segments[idx].start
        max_end = self.segments[idx].end
        
        while idx < self.nsegments:

            # skip empty
            if self.segments[idx].start == self.segments[idx].end: 
                idx += 1
                continue
            # no overlap - save segment
            if self.segments[idx].start >= max_end:
                self.segments[insertion_idx].end = max_end
                insertion_idx += 1
                assert insertion_idx < self.nsegments, "insertion_idx out of range: %i >= %i" % (insertion_idx, self.nsegments)
                self.segments[insertion_idx].start = self.segments[idx].start

            max_end = lmax( self.segments[idx].end, max_end )
            idx += 1

        # store final segment
        self.segments[insertion_idx].end = max_end

        # truncate array
        insertion_idx += 1
        self.nsegments = insertion_idx
        self._resize( self.nsegments )
        self.is_normalized = 1

    cpdef merge( self, PositionDifference distance ):
        '''merge all overlapping segments and remove empty segments.

        This function is a generalization of the merge method and
        will merge adjacent and/or neighbouring segments.

        If *distance* = 0, adjacent segments will be merged. 

        If*distance* = n, neighbouring segments with a separation of at most
        *n* bases will be merged.

        This function works in-place.
        '''

        cdef PositionDifference max_end
        cdef int idx, insertion_idx
        if self.nsegments == 0:
            self.is_normalized = 1
            return

        self.sort()

        insertion_idx = 0

        # skip to first non-empty segment
        idx = 0
        while idx < self.nsegments and self.segments[idx].start == self.segments[idx].end: 
            idx += 1

        # only empty segments - empty list, but declare normalized
        if idx == self.nsegments:
            self.nsegments = 0
            self.is_normalized = 1
            return

        self.segments[insertion_idx].start = self.segments[idx].start
        max_end = self.segments[idx].end
        
        while idx < self.nsegments:
            # skip empty
            if self.segments[idx].start == self.segments[idx].end: 
                idx += 1
                continue
            # no overlap - save segment
            # note the cange from > to >=
            if <PositionDifference>self.segments[idx].start - distance > max_end:
                self.segments[insertion_idx].end = max_end
                insertion_idx += 1
                self.segments[insertion_idx].start = self.segments[idx].start

            max_end = lmax( self.segments[idx].end, max_end )
            idx += 1

        # store final segment
        self.segments[insertion_idx].end = max_end

        # truncate array
        insertion_idx += 1
        self.nsegments = insertion_idx
        self._resize( self.nsegments )
        self.is_normalized = 1

    cpdef check( self ):
        '''check if segment list is normalized.

        If it is, return True and set flag is_normalized.
        '''

        # empty lists are normalized
        if self.nsegments == 0:
            self.is_normalized = 1
            return self.is_normalized

        self.is_normalized = 0

        cdef int idx
        cdef Segment segment

        segment = self.segments[0]
        if segment.start >= segment.end:
            raise ValueError( "empty/invalid segment in segmentlist: %s" % str( segment ) )

        for idx from 1 <= idx < self.nsegments:
            segment = self.segments[idx]
            if segment.start >= segment.end:
                raise ValueError( "empty/invalid segment in segmentlist: %s" % str( segment ) )
            if self.segments[idx-1].start > segment.start:
                raise ValueError( "segment list is not sorted: %s > %s" % (self.segments[idx-1], segment) )
            if self.segments[idx-1].end > segment.start:
                raise ValueError( "segment overlap: %s overlaps %s" % (self.segments[idx-1], segment) )
            #if self.segments[idx-1].end == segment.start:
            #    raise ValueError( "segment touch: %s touches %s" % (self.segments[idx-1], segment) )

        self.is_normalized = 1

        return self.is_normalized

    cdef int _getInsertionPoint( self, Segment other ):
        '''return insertion point for other.

        The insertion point denotes the element at which
        or after which *other* should be inserted to keep
        the sort order.

        The function returns -1, nsegments if other
        is before the first or after the last element, respectively.

        '''
        cdef int idx
        assert self.is_normalized, "searching in non-normalized list"
        if self.nsegments == 0: return -1

        # avoid out of range searches
        if other.start >= self.segments[self.nsegments-1].end:
            return self.nsegments
        if other.end <= self.segments[0].start: 
            return -1

        idx = searchsorted( self.segments,
                            self.nsegments,
                            sizeof( Segment ),
                            &other,
                            &cmpSegments )

        if idx == self.nsegments:
            return idx-1
        elif self.segments[idx].start != other.start:
            return idx-1
        else:
            return idx

    def getInsertionPoint( self, start, end ):
        return self._getInsertionPoint( Segment( start, end ) )

    property isNormalized:
        def __get__(self): return self.is_normalized

    property isEmpty:
        def __get__(self): return self.nsegments == 0

    cpdef Position getRandomPosition( self ):
        '''return a random position within the workspace.

        Not efficient, see the specialized samplers instead if
        you want to sample repeatedly from the same workspace.
        '''
        cdef Position pos = numpy.random.randint( 0, self.sum() )
        cdef idx
        cdef Position l
        for idx from 0 <= idx < self.nsegments:
            l = self.segments[idx].end - self.segments[idx].start
            if pos > l:
                pos -= l
            else:
                return self.segments[idx].start + pos
        assert False

    cdef Position overlap( self, Segment other ):
        '''return the size of intersection between
           segment list and Segment other'''

        cdef int idx
        idx = self._getInsertionPoint( other )

        # deal with border cases
        # permit partial overlap
        if idx == self.nsegments: idx -=1
        elif idx == -1: idx=0

        cdef Position count
        count = 0

        while idx < self.nsegments and self.segments[idx].start <= other.end:
            count += segment_overlap( self.segments[idx], other )
            idx += 1
        return count

    cdef SegmentList getOverlappingSegments( self, Segment other ):
        '''return the segments overlapping
           segment list and Segment other'''

        cdef int idx
        idx = self._getInsertionPoint( other )

        # deal with border cases
        # permit partial overlap
        if idx == self.nsegments: idx -=1
        elif idx == -1: idx=0

        count = 0
        cdef SegmentList result = SegmentList()
        
        while idx < self.nsegments and self.segments[idx].start <= other.end:
            result._add( self.segments[idx] )
            idx += 1
        result.normalize()
        return result

    cpdef Position overlapWithRange( self, Position start, Position end ):
        '''return the size of intersection between
           segment list and Segment other'''

        cdef Segment s
        s = Segment( start, end )
        return self.overlap( s )

    cpdef SegmentList getOverlappingSegmentsWithRange( self, Position start, Position end ):
        '''return a list with segments overlapping range.'''
        cdef Segment s
        s = Segment( start, end )
        return self.getOverlappingSegments( s )

    cpdef Position overlapWithSegments( self, SegmentList other ):
        '''return the number of nucleotides overlapping between this and *other*.'''

        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self.sum()

        cdef int this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment = Segment(0,0)
        cdef Segment other_segment = Segment(0,0)

        cdef Position overlap
        overlap = 0

        while this_idx < self.nsegments and other_idx < other.nsegments:

            # print "this=", this_idx, self.nsegments, "other=", other_idx, other.nsegments
            if last_this_idx != this_idx:
                this_segment = self.segments[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # print this_segment, other_segment
            # skip segments in this not overlapping other
            if this_segment.end <= other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_segment.start:
                other_idx += 1
            else:
                # deal with overlap
                overlap += segment_overlap_raw( this_segment, other_segment )
                if this_segment.end < other_segment.end:
                    this_idx += 1
                elif other_segment.end < this_segment.end:
                    other_idx += 1
                else:
                    this_idx += 1
                    other_idx += 1

        return overlap

    cpdef Position intersectionWithSegments( self, 
                                             SegmentList other, 
                                             mode = "base" ):
        '''return the number of segments overlapping with *other*.

        *mode* can be either ``base`` or ``midpoint``. With ``base``, an overlap is 
        recorded if a single ``base`` in self overlaps a interval in other. With ``midpoint``,
        an overlap is counted only if the midpoint of an interval overlaps any interval in ``other``.
        '''

        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self.sum()

        cdef int this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment = Segment(0,0)
        cdef Segment other_segment = Segment(0,0)

        cdef int midpoint_overlap = mode == "midpoint"

        cdef Position noverlap
        noverlap = 0

        while this_idx < self.nsegments and other_idx < other.nsegments:

            # print "this=", this_idx, self.nsegments, "other=", other_idx, other.nsegments
            if last_this_idx != this_idx:
                this_segment = self.segments[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # skip segments in this not overlapping other
            if this_segment.end <= other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_segment.start:
                other_idx += 1
            else:
                # deal with overlap
                if midpoint_overlap:
                    if other_segment.start <= this_segment.start + (this_segment.end - this_segment.start) // 2 \
                            < other_segment.end:
                        noverlap += 1
                else:
                    noverlap += 1
                this_idx += 1

        return noverlap

    def getLengthDistribution( self, bucket_size, nbuckets ):
        '''build histogram of segments lengths.'''

        cdef int idx, i
        cdef Position l
        cdef Segment s
        cdef numpy.ndarray[DTYPE_INT_t, ndim=1] histogram
        histogram = numpy.zeros( nbuckets, dtype=numpy.int )
        for idx from 0 <= idx < self.nsegments:
            l = segment_length(self.segments[idx])
            i = <int>((l+bucket_size-1)/bucket_size)
            if i >= nbuckets:
                raise ValueError( "segment %i-%i too large: %i >= %i, increase nbuckets or bucket_size such that nbuckets * bucket_size > %i" %\
                                      ( self.segments[idx].start,
                                        self.segments[idx].end,
                                        l, 
                                        nbuckets * bucket_size,
                                        l
                                        ) )
            histogram[i] += 1
        return histogram

    cdef SegmentList truncate( self, Segment other ):
        '''truncate Segment list to range given by segment.'''
        
        cdef int idx
        cdef Segment * s
        cdef Position start = other.start
        cdef Position end = other.end

        for idx from 0 <= idx < self.nsegments:
            s = & self.segments[idx] 
            if s.end < start: s.start = s.end = 0
            elif s.start > end: s.start = s.end = 0
            else:
                if s.start < start: s.start = start
                if s.end > end: s.end = end

        self.normalize()
        return self

    cpdef SegmentList getFilledSegmentsFromStart( self, Position start, PositionDifference remainder ):
        '''start filling segment from *start* until *remainder*
           bases have been covered. 

           This method wraps around at the ends if remainder can not
           be filled. If *remainder* is larger than the segment list,
           the segment list itself will be returned.

           return a segment list 
           '''

        if remainder > self.sum():
            return SegmentList( clone = self )

        cdef SegmentList result = SegmentList()

        cdef int idx
        idx = self.getInsertionPoint( start, start + 1 )
        cdef Position end
        
        # deal with border cases
        # permit partial overlap
        if idx == self.nsegments: idx -=1
        elif idx == -1: idx=0

        # add from start until complete
        while remainder > 0:
            if self.segments[idx].end < start: 
                pass
            else:
                start = lmax( self.segments[idx].start, start )
                end = lmin( self.segments[idx].end, start + remainder )
                remainder -= end - start
                result._add( Segment( start, end ) )
                
            idx += 1
            if idx == self.nsegments:
                idx = 0
                start = self.segments[idx].start

        result.normalize()
        return result

    cpdef SegmentList getFilledSegmentsFromEnd( self, Position end, PositionDifference remainder ):
        '''start filling segment from *end* until *remainder*
           bases have been covered. Fill in reverse order.

           This method wraps around at the ends if remainder can not
           be filled. If *remainder* is larger than the segment list,
           the segment list itself will be returned.

           return a segment list 
           '''

        if remainder > self.sum():
            return SegmentList( clone = self )

        cdef SegmentList result = SegmentList()

        cdef int idx
        idx = self.getInsertionPoint( end, end + 1 )
        cdef Position start
        
        # deal with border cases
        # permit partial overlap
        if idx == self.nsegments: idx -=1
        elif idx == -1: idx=0

        # add from start until complete
        while remainder > 0:
            if self.segments[idx].start > end: 
                pass
            else:
                end = lmin( self.segments[idx].end, end )
                start = lmax( self.segments[idx].start, end - remainder )
                remainder -= end - start
                result._add( Segment( start, end ) )
                
            idx -= 1
            if idx < 0:
                idx = self.nsegments - 1
                end = self.segments[idx].end

        result.normalize()
        return result

    cpdef SegmentList filter( self, SegmentList other ):
        '''remove all segments that are not in *other*

        Both this and other are assumed to have been normalized.
        '''

        # avoid self-self comparison
        if other.segments == self.segments: return self
        if self.nsegments == 0: return self

        cdef Segment * new_segments
        cdef size_t allocated
        # create new list
        new_segments =<Segment*>malloc( self.nsegments * sizeof( Segment ) )

        cdef int this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        working_idx = this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment = Segment(0,0)
        cdef Segment other_segment = Segment(0,0)
        # for negative segments, do not use -1
        cdef Position last_start = self.segments[0].start - 1

        while this_idx < self.nsegments and other_idx < other.nsegments:

            # print "this=", this_idx, self.nsegments, "other=", other_idx, other.nsegments
            if last_this_idx != this_idx:
                this_segment = self.segments[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # print this_segment, other_segment
            # skip segments in this not overlapping other
            if this_segment.end <= other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_segment.start:
                other_idx += 1
            else:
                # deal with overlap
                if last_start != this_segment.start:
                    new_segments[working_idx].start = this_segment.start
                    new_segments[working_idx].end = this_segment.end
                    working_idx += 1
                    last_start = this_segment.start

                if this_segment.end < other_segment.end:
                    this_idx += 1
                elif other_segment.end < this_segment.end:
                    other_idx += 1
                else:
                    this_idx += 1
                    other_idx += 1

        free( self.segments )
        self.nsegments = working_idx
        self.segments = new_segments
        self._resize( self.nsegments )
        return self

    cpdef SegmentList intersect( self, SegmentList other ):
        '''intersect this segment list with another.

        Intervals are truncated to the overlap between this and
        other.

        Both this and other are assumed to have been normalized.

        Segments are truncated to the smallest intersection. Thus,
        the intersection between ``[(0,10), (10,20)]`` and ``[(0,20)]``
        will be ``[(0,10), (10,20)]`` and not ``[(0,20)]``.

        '''
        assert self.is_normalized, "intersection of a non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self

        if self.nsegments == 0: return self

        cdef Segment * new_segments
        cdef size_t allocated
        # create new list with 10% overhead
        allocated = int(lmax( self.nsegments, other.nsegments) * 1.1)
        new_segments =<Segment*>malloc( allocated * sizeof( Segment ) )

        cdef int this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        working_idx = this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment = Segment(0,0)
        cdef Segment other_segment = Segment(0,0)

        while this_idx < self.nsegments and other_idx < other.nsegments:

            # print "this=", this_idx, self.nsegments, "other=", other_idx, other.nsegments
            if last_this_idx != this_idx:
                this_segment = self.segments[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # print this_segment, other_segment
            # skip segments in this not overlapping other
            if this_segment.end <= other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_segment.start:
                other_idx += 1
            else:
                # deal with overlap
                new_segments[working_idx].start = lmax(this_segment.start, other_segment.start )
                new_segments[working_idx].end = lmin(this_segment.end, other_segment.end )
                working_idx += 1
                if working_idx >= allocated:
                    allocated *= 2
                    new_segments = <Segment*>realloc( new_segments, allocated * sizeof(Segment ) )
                    assert new_segments != NULL

                if this_segment.end < other_segment.end:
                    this_idx += 1
                elif other_segment.end < this_segment.end:
                    other_idx += 1
                else:
                    this_idx += 1
                    other_idx += 1

        free( self.segments )
        self.segments = new_segments
        self.nsegments = working_idx
        self.allocated = allocated
        return self

    cpdef extend_segments( self, int extension ):
        '''extend all intervals by a certain amount.

        If the extension causes an interval to extend beyond the chromosome
        start, the interval is truncated.

        The segment list will not be normalized afterwards.'''

        cdef int idx
        cdef Segment * s
        for idx from 0 <= idx < self.nsegments:
            s = &self.segments[idx]
            s.start -= lmin( extension, s.start)
            s.end += extension
        self.is_normalized = False

    cpdef expand_segments( self, double expansion ):
        '''expand all intervals by a certain amount.

        The interval length is multiplied by *expansion* and 
        extended around the mid-point.

        If the expansion causes an interval to extend beyond the chromosome
        start, the interval is truncated.

        The segment list will not be normalized afterwards.'''

        cdef int idx, 
        cdef PositionDifference l, extension
        cdef Segment * s
        cdef double e_1 = (expansion - 1.0) / 2.0
        
        if expansion <= 0: raise ValueError( "invalid expansion: %f <= 0" % expansion)

        for idx from 0 <= idx < self.nsegments:
            s = &self.segments[idx]
            l = s.end - s.start 
            extension = <PositionDifference>floor(l * e_1 )
            s.start -= lmin( extension, s.start)
            s.end += extension
        self.is_normalized = False

    def shift( self, PositionDifference offset ):
        '''shift segments by a certain offset.
        
        raises ValueError if segment coordinates become negative.
        '''
        cdef int idx
        cdef Segment * s
        for idx from 0 <= idx < self.nsegments:
            s = &self.segments[idx]
            s.end += offset
            s.start += offset
            if s.start < 0: raise ValueError( "shift creates negative coordinates" )

    cpdef Position sum( self ):
        '''return total length of all segments.'''
        cdef Position total
        cdef int idx
        cdef Segment s
        total = 0
        for idx from 0 <= idx < self.nsegments:
            s = self.segments[idx]
            total += s.end - s.start
        return total

    cpdef Position max( self ):
        '''return maximum coordinate.'''
        assert self.is_normalized, "maximum from non-normalized list"
        if self.nsegments == 0: return 0
        return self.segments[self.nsegments - 1].end

    cpdef Position min( self ):
        '''return minimum coordinate.'''
        assert self.is_normalized, "minimum from non-normalized list"
        if self.nsegments == 0: return 0
        return self.segments[0].start

    cpdef clear( self ):
        '''removes all segments from list.

        Note that this method does not shrink the allocated memory.
        '''
        self.nsegments = 0
        self.is_normalized = 1

    def __len__(self):
        return self.nsegments

    def __dealloc__(self):
        if self.segments != NULL:
            free( self.segments )

    def __str__(self):
        return str(self.asList())

    def asList( self ):

        cdef int idx
        result = []
        for idx from 0 <= idx < self.nsegments:
            result.append( (self.segments[idx].start, self.segments[idx].end) )
        return result

    def asLengths( self ):
        cdef int idx
        result = []
        for idx from 0 <= idx < self.nsegments:
            result.append( self.segments[idx].end - self.segments[idx].start)
        return result

    def __iter__( self ):
        return SegmentListIterator( self )

    def __getitem__(self, key ):
        if key >= self.nsegments:
            raise IndexError( "index out of range" )
        return self.segments[key].start, self.segments[key].end

    def __cmp__(self, SegmentList other ):
        cdef int idx
        x = self.__len__().__cmp__(len(other))
        if x != 0: return x
        for idx from 0 <= idx < self.nsegments:
            x = cmpSegmentsStartAndEnd( &self.segments[idx], &other.segments[idx] )
            if x != 0: return x
        return 0

# SegmentTuple = collections.namedtuple( "SegmentTuple", "start, end" )

cdef class SegmentListIterator:

    cdef SegmentList segment_list
    cdef Position idx

    def __init__(self, SegmentList container ):
        self.idx = 0
        self.segment_list = container

    def __iter__(self): return self

    def __next__(self):
        cdef Position t
        cdef Segment v
        if self.idx >= self.segment_list.nsegments:
            raise StopIteration
        v = self.segment_list.segments[self.idx]
        self.idx += 1
        return v.start, v.end

