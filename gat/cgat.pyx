# cython: embedsignature=True
# cython: profile=True

import types, collections


cimport cython

cdef extern from "string.h":
    ctypedef int size_t
    void *memcpy(void *dest, void *src, size_t n)
    char *strtok_r(char *str, char *delim, char **saveptr)
    char *strncpy(char *dest, char *src, size_t n)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    int c_abs "abs" (int)
    int atoi( char *nptr)
    long atol( char *nptr)
    double atof( char *nptr)
    void qsort(void *base, size_t nmemb, size_t size,
               int(*compar)(void *, void *))
    int rand()
    int rand_r(unsigned int *seedp)
    void srand(unsigned int seed)

cdef extern from "stdint.h":
    ctypedef int int64_t
    ctypedef int int32_t
    ctypedef int uint32_t
    ctypedef int uint8_t
    ctypedef int uint64_t

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    char *fgets(char *str, int size, FILE *ifile)
    void fprintf(FILE * f, char* s, char* s)
    void printf(char* s,...)
    int feof(FILE *stream)
    size_t strlen(char *s)
    size_t getline(char **lineptr, size_t *n, FILE *stream)
    char *strstr(char *, char *)
    char *strchr(char *string, int c)
    int fileno(FILE *stream)

# Next, enter the builtin file class into the namespace:
#cdef extern from "fileobject.h":
#    ctypedef class __builtin__.file [object PyFileObject]:
#        pass

cdef extern from "gat_utils.h":
    long searchsorted(void * base,
                       size_t nmemb,
                       size_t size,
                       void * target,
                       int(*compar)(void *, void *))


#####################################################
## numpy import
# note both import and cimport are necessary
import numpy 
cimport numpy
DTYPE_INT = numpy.int 
ctypedef numpy.int_t DTYPE_INT_t
#DTYPE_FLOAT = numpy.float 
#ctypedef numpy.float_t DTYPE_FLOAT_t

cdef struct Segment:
    long start
    long end

# min/max are not optimized, so declare them as C functions
@cython.profile(False)
cdef inline long lmin( long a, long b):
    if a < b: return a
    return b

# min/max are not optimized, so declare them as C functions
@cython.profile(False)
cdef inline long lmax( long a, long b):
    if a > b: return a
    return b

@cython.profile(False)
cdef inline long segment_overlap(  Segment a, Segment b ):
    return lmax(0, lmin( a.end, b.end) - lmax(a.start, b.start))

@cython.profile(False)
cdef inline long segment_overlap_raw(  Segment a, Segment b ):
    return lmin(a.end, b.end) - lmax(a.start, b.start)

@cython.profile(False)
cdef inline long segment_length(  Segment a ):
    return a.end - a.start

@cython.profile(False)
cdef int cmpSegments( void * s1, void * s2 ):
    return (<Segment *>s1).start - (<Segment *>s2).start;

@cython.profile(False)
cdef int cmpLong( void * s1, void * s2 ):
    return (<long*>s1)[0] - (<long*>s2)[0]

cdef class SegmentList:
    '''list of segments.
    
    A list of (non-overlapping) segments. 

    The intersection algorithms assume that segments 
    are non-overlapping and sorted by start.

    Call normalize to merge overlapping segments.
    '''
    cdef Segment * segments
    cdef size_t nsegments
    cdef size_t allocated
    cdef int is_normalized 

    def __init__(self, long allocate = 1000, 
                 SegmentList clone = None,
                 iter = None,
                 normalize = False ):
        '''create empty list of segments.'''
        cdef long idx
        self.nsegments = 0
        # an empty list is normalized
        self.is_normalized = 1

        if clone != None:
            self.nsegments = self.allocated = clone.nsegments
            self.segments = <Segment*>calloc( clone.nsegments, sizeof( Segment ) )
            memcpy( self.segments, 
                    clone.segments,
                    clone.nsegments * sizeof(Segment ) )
            self.is_normalized = clone.is_normalized
        elif iter:
            a = tuple(iter)
            self.nsegments = self.allocated = len(a)
            self.segments = <Segment*>calloc( self.nsegments, sizeof( Segment ) )
            idx = 0
            for start, end in a:
                self.segments[idx] = Segment( start, end )
                idx += 1
            self.is_normalized = 0
        else:
            self.allocated = allocate
            self.segments = <Segment*>calloc( allocate, sizeof( Segment ) )

        if normalize: self.normalize()

    def sort( self ):
        '''sort segments.'''
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
        if new_size >= self.allocated:
            self.segments = <Segment*>realloc( <void *>self.segments, new_size * sizeof( Segment ) )
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
        
        if self.nsegments == self.allocated:
            self.allocated *= 2
            self.segments = <Segment*>realloc( self.segments, 
                                               self.allocated * sizeof( Segment ) )
        self.segments[self.nsegments] = segment
        self.nsegments += 1
        self.is_normalized = 0

    cpdef add( self, long start, long end ):
        cdef Segment segment
        segment = Segment( start, end)
        self._add( segment )

    cpdef trim( self, long pos, long size ):
        '''trim segment list by removing *size* nucleotides from
        the segment that includes *pos*.
        '''
        assert self.is_normalized, "trimming in non-normalized list"

        cdef long idx
        cdef Segment other, seg

        assert self.sum() > size, "trimming more than the total length (%i < %i)" % (self.sum(), size)
        #print "at start=%i, removing %i" % (self.sum(), size)

        other = Segment( pos, pos + 1)
        idx = searchsorted( self.segments, 
                            self.nsegments,
                            sizeof( Segment ),
                            &other,
                            &cmpSegments )
        
        # wrap around
        if idx == self.nsegments: idx = 0
            
        while size > 0:
            seg = self.segments[idx]
            if segment_length( seg ) < size:
                self.segments[idx] = Segment(0,0)
                size -= len(seg)
            else:
                self.segments[idx] = Segment( seg.start + size, seg.end )
                size = 0

            idx += 1
            if idx == self.nsegments: idx = 0

        #print "at end=%i, removing %i" % (self.sum(), size)

    cpdef normalize( self ):
        '''merge all overlapping segments.
        
        Adjacent segments are not merged.

        This function works in-place.
        '''

        cdef long idx, max_end, insertion_idx

        self.sort()
 
        insertion_idx = 0
        idx = 0
        max_end = self.segments[idx].end
        
        for idx from 0 <= idx < self.nsegments:
            if self.segments[idx].start >= max_end:
                self.segments[insertion_idx].end = max_end
                insertion_idx += 1
                self.segments[insertion_idx].start = self.segments[idx].start
            max_end = lmax( self.segments[idx].end, max_end )

        # store final segment
        self.segments[insertion_idx].end = max_end

        # truncate array
        idx += 1
        insertion_idx += 1
        self.nsegments = insertion_idx
        self.segments = <Segment*>realloc( self.segments, self.nsegments * sizeof( Segment ) )
        self.allocated = self.nsegments
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

        cdef long idx
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

    cdef long getInsertionPoint( self, Segment other ):
        '''return insertion point for other.

        The insertion point denotes the element at which
        or after which *other* should be inserted to keep 
        the sort order.
        '''
        cdef long idx
        assert self.is_normalized, "searching in non-normalized list"

        # avoid out of range searches
        if other.start > self.segments[self.nsegments-1].end:
            return self.nsegments

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

    property isNormalized:
        def __get__(self): return self.is_normalized

    property isEmpty:
        def __get__(self): return self.nsegments == 0


    cdef long overlap( self, Segment other ):
        '''return the size of intersection between 
           segment list and Segment other'''
        
        cdef long idx
        idx = self.getInsertionPoint( other )
        if idx == self.nsegments: return 0
        
        cdef long count
        count = 0

        while idx < self.nsegments and self.segments[idx].start <= other.end:
            count += segment_overlap( self.segments[idx], other )
            idx += 1
        return count

    cdef long _overlapWithRange( self, long start, long end ):
        '''return the size of intersection between 
           segment list and Segment other'''
        
        cdef Segment s
        s = Segment( start, end )
        return self.overlap( s )

    def overlapWithSegments( self, SegmentList other ):
        '''return the number of nucleotides overlapping between this and *other*.'''
        
        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self.sum()

        cdef long this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment, other_segment
        cdef long overlap
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

    def overlapWithRange( self, start, end ):
        return self._overlapWithRange( start, end )

    def getLengthDistribution( self, bucket_size, nbuckets ):
        '''build histogram of segments lengths.'''
        
        cdef long idx
        cdef long l 
        cdef Segment s
        cdef numpy.ndarray[DTYPE_INT_t, ndim=1] histogram
        histogram = numpy.zeros( nbuckets, dtype=numpy.int )
        for idx from 0 <= idx < self.nsegments:
            l = segment_length(self.segments[idx])
            histogram[(int)((l+bucket_size-1)/bucket_size)] += 1
        return histogram

    cpdef SegmentList intersect( self, SegmentList other ):
        '''intersect this segment list with another.
        
        Intervals are truncated to the overlap between this and 
        other.

        Both this and other are assumed to have been normalized.
        
        Segments are truncated to the smallest intersection. Thus,
        the intersection between ``[(0,10), (10,20)]`` and ``[(0,20)]`` 
        will be ``[(0,10), (10,20)]`` and not ``[(0,20)]``.

        '''
        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self

        cdef Segment * new_segments
        cdef size_t allocated
        # create new list with 10% overhead
        allocated = int(lmax( self.nsegments, other.nsegments) * 1.1)
        new_segments =<Segment*>malloc( allocated * sizeof( Segment ) )

        cdef long this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        working_idx = this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment, other_segment

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

    cpdef long sum( self ): 
        '''return total length of all segments.'''
        cdef long total
        cdef long idx
        cdef Segment s
        total = 0
        for idx from 0 <= idx < self.nsegments:
            s = self.segments[idx]
            total += s.end - s.start
        return total

    cpdef clear( self ):
        '''removes all segments from list.
        
        Note that this method does not shrink the allocated memory.
        '''
        self.nsegments = 0
        self.is_normalized = 0

    def __len__(self):
        return self.nsegments

    def __dealloc__(self):
        if self.segments != NULL:
            free( self.segments )

    def __str__(self):
        return str(self.asList())

    def asList( self ):

        cdef long idx
        result = []
        for idx from 0 <= idx < self.nsegments:
            result.append( (self.segments[idx].start, self.segments[idx].end) )
        return result

    def __iter__( self ):
        return SegmentListIterator( self )

    def __getitem__(self, key ):
        if key >= self.nsegments:
            raise IndexError( "index out of range" )
        return self.segments[key].start, self.segments[key].end

# SegmentTuple = collections.namedtuple( "SegmentTuple", "start, end" )

cdef class SegmentListIterator:

    cdef SegmentList segment_list
    cdef long idx 

    def __init__(self, SegmentList container ):
        self.idx = 0
        self.segment_list = container
        
    def __iter__(self): return self

    def __next__(self):
        cdef long t
        cdef Segment v
        if self.idx >= self.segment_list.nsegments:
            raise StopIteration
        v = self.segment_list.segments[self.idx]
        self.idx += 1
        return v.start, v.end

cdef class SegmentListSamplerSlow:
   
    cdef SegmentList segment_list
    cdef numpy.ndarray  cdf
    cdef long total_size

    def __init__(self, SegmentList segment_list ):
        
        assert len(segment_list) > 0, "sampling from empty segment list" 
        
        self.segment_list = segment_list
        self.cdf = numpy.cumsum( [x[1] - x[0] for x in self.segment_list ] )
        self.total_size = self.cdf[-1]

    cpdef sample( self, long length ):
        '''return a new position within segment list.'''
        
        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef long r, offset, pos, overlap
        cdef size_t segment_index

        r = numpy.random.randint( 0, self.total_size )
        segment_index = numpy.searchsorted( self.cdf, r )
        offset = r - self.cdf[segment_index]
        if offset == 0:
            pos = self.segment_list.segments[segment_index].start
        else:
            pos = self.segment_list.segments[segment_index].end + offset

        overlap = self.segment_list._overlapWithRange( pos, pos+length )

        assert overlap > 0, "sample %i does not overlap with workspace: offset=%i, r=%i, index=%i/%i" %\
            (pos, r - self.cdf[segment_index], r, segment_index, self.segment_list.nsegments)

        return pos, overlap

cdef class SegmentListSampler:
   
    cdef SegmentList segment_list
    cdef long * cdf
    cdef long total_size
    cdef long nsegments

    def __init__(self, SegmentList segment_list ):
        cdef long i, totsize
        assert len(segment_list) > 0, "sampling from empty segment list" 
        
        self.segment_list = segment_list
        self.nsegments = len(segment_list)
        self.cdf = <long*>malloc( sizeof(long) * self.nsegments )
        self.total_size = 0
        for i from 0 <= i < len(segment_list):
            self.total_size += segment_length( segment_list.segments[i] )
            self.cdf[i] = self.total_size

    cpdef sample( self, long length ):
        '''return a new position within segment list.'''
        
        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef long r, offset, pos, overlap
        cdef size_t segment_index

        # r = rand() / (RAND_MAX / total_size + 1)

        r = numpy.random.randint( 0, self.total_size )
        segment_index = searchsorted( self.cdf, 
                                      self.nsegments,
                                      sizeof(long),
                                      &r,
                                      &cmpLong,
                                      )

        offset = r - self.cdf[segment_index]
        if offset == 0:
            pos = self.segment_list.segments[segment_index].start
        else:
            pos = self.segment_list.segments[segment_index].end + offset

        overlap = self.segment_list._overlapWithRange( pos, pos+length )

        assert overlap > 0, "sample %i does not overlap with workspace: offset=%i, r=%i, index=%i/%i" %\
            (pos, r - self.cdf[segment_index], r, segment_index, self.segment_list.nsegments)

        return pos, overlap

cdef class HistogramSamplerSlow:

    cdef numpy.ndarray cdf
    cdef long bucket_size

    def __init__(self, 
                 numpy.ndarray histogram,
                 long bucket_size ):
        
        assert len(histogram) > 0, "sampling from empty histogram" 

        self.cdf = numpy.cumsum( histogram, dtype = numpy.float )
        self.cdf /= self.cdf[-1]
        self.bucket_size = bucket_size

    cpdef long sample( self ):
        '''return a new position within segment list.'''
        
        cdef long base, ip
        cdef double r
        cdef long bucket_size
        bucket_size = self.bucket_size
        # note: could be made quicker by
        # - creating a sample random doubles once, see numpy.random_integers?
        # - avoiding the binary search?
        r = numpy.random.random_sample( 1 )[0]
        ip = numpy.searchsorted( self.cdf, r )
        base = ip*bucket_size
        if bucket_size>1 : return base + numpy.randint(0, bucket_size)

        return base

cdef class HistogramSampler:

    cdef long * cdf
    cdef long bucket_size
    cdef long nbuckets
    cdef long total_size

    @cython.boundscheck(False)
    def __init__(self, 
                 numpy.ndarray[DTYPE_INT_t, ndim=1] histogram,
                 long bucket_size ):
        cdef long i
        self.cdf = NULL

        self.nbuckets = len(histogram)
        assert self.nbuckets > 0, "sampling from empty histogram" 

        self.cdf = <long*>malloc( sizeof(long) * self.nbuckets )
        self.total_size = 0
        for i from 0 <= i < self.nbuckets:
            self.total_size += histogram[i]
            self.cdf[i] = self.total_size
        
        self.bucket_size = bucket_size

    cpdef long sample( self ):
        '''return a new position within segment list.'''
        
        cdef long base, index, r
        
        r = numpy.random.randint( 0, self.total_size )
        index = searchsorted( self.cdf, 
                              self.nbuckets,
                              sizeof(long),
                              &r,
                              &cmpLong,
                              )

        # + 1 to avoid 0 length
        base = index * self.bucket_size + 1

        if self.bucket_size > 1 : 
            return base + numpy.randint(0, self.bucket_size)

        return base

    def __dealloc__(self ):
        if self.cdf != NULL:
            free( self.cdf )
cdef class Sampler:
    pass

cdef class SamplerAnnotator(Sampler):

    cdef long bucket_size
    cdef long nbuckets
    cdef long nunsuccessful_rounds

    def __init__( self, bucket_size = 1, nbuckets = 100000 ):
        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.nunsuccessful_rounds

    cpdef SegmentList sample( self,
                              SegmentList segments, 
                              SegmentList workspace ):
        '''main sampling procedure.'''

        cdef long remaining
        cdef long true_remaining
        cdef int nunsuccessful_rounds 
        cdef int max_unsuccessful_rounds
        cdef long length
        cdef long pos
        cdef long overlap
        cdef Segment segment, sampled_segment
        cdef SegmentList sampled_segments
        cdef SegmentList unintersected_segments
        cdef SegmentList intersected_segments
        cdef SegmentList working_segments
        cdef SegmentListSampler sls, temp_sampler

        assert segments.is_normalized, "segment list is not normalized"
        assert workspace.is_normalized, "workspace is not normalized"

        unintersected_segments = SegmentList()
        intersected_segments = SegmentList()

        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )
        working_segments.intersect( workspace )
        
        if len(working_segments) == 0:
            return intersected_segments

        # create space for sampled segments, add additional 10%
        # safety margin to avoid realloc calls
        sampled_segments = SegmentList( allocate = int(1.1 * len(segments)) )

        # build length histogram
        histogram = working_segments.getLengthDistribution( self.bucket_size, 
                                                            self.nbuckets )

        hs = HistogramSampler( histogram, self.bucket_size )

        # get segment sampler
        sls = SegmentListSampler( workspace )

        ltotal = working_segments.sum()
        remaining = ltotal
        true_remaining = remaining
        nunsuccessful_rounds = 0
        max_unsuccessful_rounds = 20

        while true_remaining > 0 and nunsuccessful_rounds < max_unsuccessful_rounds:

            # print "------------------------------------"
            # print "true_remanining", true_remaining, "remaining", remaining, "ltotal", ltotal
            # print "sampled_segments", sampled_segments.sum(), len(sampled_segments)
            # print "unintersected_segments", unintersected_segments.sum(), len(unintersected_segments)
            # print "intersected_segments", intersected_segments.sum(), len(intersected_segments)

            # Sample a segment length from the histogram
            length = hs.sample()
            assert length > 0

            # If we can now potentially get required number of nucleotides, recompute
            # the required amount since we might need more due to overlap
            if remaining <= length:
                # print "recomputing segments"
                # print " sampled segments", sampled_segments.asList()
                unintersected_segments.extend( sampled_segments )
                # print "extended"
                sampled_segments.clear()
                # compute overlap with workspace
                intersected_segments = SegmentList( clone = unintersected_segments )
                intersected_segments.normalize()
                intersected_segments.intersect( workspace )
                # print " sampled_segments", sampled_segments.sum(), len(sampled_segments)
                # print " unintersected_segments", unintersected_segments.sum(), len(unintersected_segments)
                # print unintersected_segments.asList()
                # print " intersected_segments", intersected_segments.sum(), len(intersected_segments)
                # print intersected_segments.asList()
                # assert intersected_segments.sum() <= unintersected_segments.sum()
                # if not intersected_segments.sum() <= ltotal:
                #     print "\n".join( [str(x) for x in intersected_segments ] ) + "\n"
                # assert intersected_segments.sum() <= ltotal, \
                #     "length of segments exceeds ltotal: %i > %i" % (intersected_segments.sum(), ltotal)
                remaining = ltotal - intersected_segments.sum()
                if true_remaining == remaining:
                    nunsuccessful_rounds += 1
                else:
                    true_remaining = remaining

            # deal with overshoot
            if true_remaining < 0:
                # print "removing overshoot"
                temp_sampler = SegmentListSampler( unintersected_segments )
                pos, overlap = temp_sampler.sample( length )
                unintersected_segments.trim( pos, -true_remaining )

            # sample a position until we get a nonzero overlap
            pos, overlap = sls.sample( length )
            sampled_segment = Segment( pos, pos+length )

            # If intersect > remaining we could well add too much (but we're not sure, since
            # there may be overlap with previously sampled segments).  Make sure we don't
            # overshoot the mark.  However, if we can't increase the amount of overlap in a
            # reasonable amount of time, don't bother
            if remaining < overlap and nunsuccessful_rounds < max_unsuccessful_rounds:
                # print "truncating sampled segment", remaining, overlap
                if workspace._overlapWithRange( pos, pos+1 ) > 0:
                    # Left end overlaps with workspace
                    sampled_segment.end = sampled_segment.start + remaining;
                else:
                    # Right end does?
                    sampled_segment.start = sampled_segment.end - remaining;
                overlap = workspace._overlapWithRange( sampled_segment.start, sampled_segment.end );

            if true_remaining > 0:
                # print "adding segment ", sampled_segment, true_remaining
                sampled_segments._add( sampled_segment )
                remaining -= overlap

        self.nunsuccessful_rounds = nunsuccessful_rounds
        return intersected_segments


cdef class Counter:
    '''base class for objects that compute counts
    between two segment lists.
    '''

cdef class CounterNucleotides(Counter):
    def __init__(self): 
        pass

    def __call__(self, segments, annotations, workspace = None ):
        '''return number of nucleotides overlapping between segments and annotations.'''
        return segments.overlapWithSegments( annotations )

cdef class CounterDensities(Counter):
    '''count overlap densities.

    In order to compare with Annotator results. The number of bases
    overlapping are return as densities
    '''

    def __init__(self): 
        pass

    def __call__(self, segments, annotations, workspace ):
        '''return number of nucleotides overlapping between segments and annotations.'''
        return segments.overlapWithSegments( annotations )
        
        
        
