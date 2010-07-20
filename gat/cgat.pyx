# cython: embedsignature=True
# cython: profile=True

import types, collections, re, os
import IOTools
import Experiment as E

cimport cython

cdef extern from "string.h":
    ctypedef int size_t
    void *memcpy(void *dest, void *src, size_t n)
    char *strtok_r(char *str, char *delim, char **saveptr)
    char *strncpy(char *dest, char *src, size_t n)
    void *memchr(void *s, int c, size_t n)

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

cdef extern from "stdio.h":
    ctypedef struct FILE
    # check 32/64 bit compatibility for large files
    ctypedef long off_t
    ctypedef struct fpos_t:
        pass
    cdef int SEEK_SET
    cdef int SEEK_CUR
    cdef int SEEK_END
    FILE *fopen(char *,char *)
    int fclose(FILE *)
    int feof(FILE *)
    int sscanf(char *str,char *fmt,...)
    int sprintf(char *str,char *fmt,...)
    int fprintf(FILE *ifile,char *fmt,...)
    int ferror(FILE *stream)
    int fseeko(FILE *stream, off_t offset, int whence)
    off_t ftello(FILE *stream)
    size_t fwrite( void *ptr,
                   size_t size,
                   size_t nmemb,
                   FILE *stream)
    size_t fread(void *ptr,
                 size_t size,
                 size_t nmemb,
                 FILE *stream)
    char *fgets(char *str,int size,FILE *ifile)
    int fgetpos(FILE *stream, fpos_t *pos)
    int fsetpos(FILE *stream, fpos_t *pos)
    int printf(char *format, ...)

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    char *fgets(char *str, int size, FILE *ifile)
    int feof(FILE *stream)
    size_t strlen(char *s)
    size_t getline(char **lineptr, size_t *n, FILE *stream)
    char *strstr(char *, char *)
    char *strchr(char *string, int c)
    int fileno(FILE *stream)

# Next, enter the builtin file class into the namespace:
cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

cdef extern from "gat_utils.h":
    long searchsorted(void * base,
                       size_t nmemb,
                       size_t size,
                       void * target,
                       int(*compar)(void *, void *))

    long searchargsorted(void * base,
                         int * sorted,
                         size_t nmemb,
                         size_t size,
                         void * target,
                         int(*compar)(void *, void *))

    int toCompressedFile( unsigned char *, size_t, FILE *)
    int fromCompressedFile( unsigned char *, size_t, FILE *)

#####################################################
## numpy import
# note both import and cimport are necessary
import numpy
cimport numpy
DTYPE_INT = numpy.int
ctypedef numpy.int_t DTYPE_INT_t
DTYPE_FLOAT = numpy.float
ctypedef numpy.float_t DTYPE_FLOAT_t

cdef struct Segment:
    long start
    long end

# min/max are not optimized, so declare them as C functions
@cython.profile(False)
cdef inline long lmin( long a, long b):
    if a < b: return a
    return b
@cython.profile(False)
cdef inline double dmin( double a, double b):
    if a < b: return a
    return b

# min/max are not optimized, so declare them as C functions
@cython.profile(False)
cdef inline long lmax( long a, long b):
    if a > b: return a
    return b
@cython.profile(False)
cdef inline double dmax( double a, double b):
    if a > b: return a
    return b

@cython.profile(False)
cdef inline long segment_overlap(  Segment a, Segment b ):
    return lmax(0, lmin( a.end, b.end) - lmax(a.start, b.start))

@cython.profile(False)
cdef inline long range_overlap( long astart, long aend, long bstart, long bend ):
    return lmax(0, lmin( aend, bend) - lmax(astart, bstart))

@cython.profile(False)
cdef inline long segment_overlap_raw(  Segment a, Segment b ):
    return lmin(a.end, b.end) - lmax(a.start, b.start)

@cython.profile(False)
cdef inline long segment_length(  Segment a ):
    return a.end - a.start

# trick to permit const void * in function definitions
cdef extern from *:
    ctypedef void * const_void_ptr "const void*"

@cython.profile(False)
cdef int cmpSegments( const_void_ptr s1, const_void_ptr s2 ):
    return (<Segment *>s1).start - (<Segment *>s2).start;

@cython.profile(False)
cdef int cmpLong( const_void_ptr s1, const_void_ptr s2 ):
    return (<long*>s1)[0] - (<long*>s2)[0]

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
        if self.nsegments == 0: return

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
        '''merge all overlapping segments and remove empty segments.

        Adjacent segments are not merged.

        This function works in-place.
        '''

        cdef long idx, max_end, insertion_idx
        if self.nsegments == 0:
            self.is_normalized = 1
            return

        self.sort()

        insertion_idx = 0
        idx = 0
        while idx < self.nsegments and self.segments[idx].start == self.segments[idx].end: 
            idx += 1
        # only empty segments - empty list, but declare normalized
        if idx == self.nsegments:
            self.nsegments = 0
            self.is_normalized = 1
            return

        max_end = self.segments[idx].end

        for idx from 0 <= idx < self.nsegments:
            if self.segments[idx].start == self.segments[idx].end: continue
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
        if self.nsegments == 0: return 0

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

        # deal with border cases
        # permit partial overlap
        if idx == self.nsegments: idx -= 1
        elif idx == -1: idx=0

        cdef long count
        count = 0

        while idx < self.nsegments and self.segments[idx].start <= other.end:
            count += segment_overlap( self.segments[idx], other )
            idx += 1
        return count

    cpdef long overlapWithRange( self, long start, long end ):
        '''return the size of intersection between
           segment list and Segment other'''

        cdef Segment s
        s = Segment( start, end )
        return self.overlap( s )

    cpdef long overlapWithSegments( self, SegmentList other ):
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

    cpdef long intersectionWithSegments( self, SegmentList other ):
        '''return the number of segments overlapping with *other*.'''

        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self.sum()

        cdef long this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment, other_segment
        cdef long noverlap
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
                noverlap += 1
                this_idx += 1

        return noverlap

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

    cpdef long max( self ):
        '''return maximum coordinate.'''
        assert self.is_normalized, "maximum from non-normalized list"
        if self.nsegments == 0: return 0
        return self.segments[self.nsegments - 1].end

    cpdef long min( self ):
        '''return minimum coordinate.'''
        assert self.is_normalized, "minimum from non-normalized list"
        if self.nsegments == 0: return 0
        return self.segments[0].start

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

        overlap = self.segment_list.overlapWithRange( pos, pos+length )

        assert overlap > 0, "sample %i does not overlap with workspace: offset=%i, r=%i, index=%i/%i" %\
            (pos, r - self.cdf[segment_index], r, segment_index, self.segment_list.nsegments)

        return pos, overlap

cdef class SegmentListSamplerWithEdgeEffects:

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
        '''return a new position within segment list.

        This method both samples a position within the workspace
        and a position within length.
        '''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef long rpos, rseg, offset, pos, overlap, start, end
        cdef size_t segment_index

        # r = rand() / (RAND_MAX / total_size + 1)
        rpos = numpy.random.randint( 0, self.total_size )
        segment_index = searchsorted( self.cdf,
                                      self.nsegments,
                                      sizeof(long),
                                      &rpos,
                                      &cmpLong,
                                      )

        offset = rpos - self.cdf[segment_index]
        if offset == 0:
            pos = self.segment_list.segments[segment_index+1].start
        else:
            pos = self.segment_list.segments[segment_index].end + offset

        rseg = numpy.random.randint( 0, length )

        start, end = pos - rseg, pos - rseg + length

        # print "sample %i-%i, offset=%i, pos=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
        #     (start, end, offset, pos, rpos, rseg,
        #      segment_index, self.segment_list.nsegments,
        #      self.segment_list.segments[segment_index] )

        assert overlap > 0, "sample %i-%i does not overlap with workspace: offset=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
            (start, end, offset, rpos, rseg,
             segment_index, self.segment_list.nsegments,
             self.segment_list.segments[segment_index])

        return start, end, overlap

cdef class SegmentListSampler:
    '''return a new position within segment list.

    The probability of overlap is proportional
    to the workspace size.

    In order to avoid edge effects the position
    of the sample within a workspace segment is
    sampled again.
    '''

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
            # -1, as searchsorted inserts on left
            self.cdf[i] = self.total_size -1

    cpdef sample( self, long length ):
        '''return a sample.'''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef long rpos, rseg, offset, pos, overlap, start, end, lseg
        cdef size_t segment_index
        cdef Segment s

        # r = rand() / (RAND_MAX / total_size + 1)
        rpos = numpy.random.randint( 0, self.total_size )

        segment_index = searchsorted( self.cdf,
                                      self.nsegments,
                                      sizeof(long),
                                      &rpos,
                                      &cmpLong,
                                      )

        offset = rpos - self.cdf[segment_index]

        s = self.segment_list.segments[segment_index]
        start = s.start - length + 1
        lseg = segment_length( s )

        rseg = numpy.random.randint( 0, lseg + length - 1)
        start += rseg
        end = start + length

        # print "sample %i-%i, offset=%i, pos=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
        #     (start, end, offset, pos, rpos, rseg,
        #      segment_index, self.segment_list.nsegments,
        #      self.segment_list.segments[segment_index] )
        overlap = range_overlap( s.start, s.end, start, end )

        assert overlap > 0, "sample %i-%i does not overlap with workspace: offset=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
            (start, end, offset, rpos, rseg,
             segment_index,
             self.segment_list.nsegments,
             self.segment_list.segments[segment_index])

        return start, end, overlap

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

        # 1 to avoid 0 length
        if self.total_size > 1:
            r = numpy.random.randint( 1, self.total_size )
        else:
            r = 1

        index = searchsorted( self.cdf,
                              self.nbuckets,
                              sizeof(long),
                              &r,
                              &cmpLong,
                              )
        base = index * self.bucket_size

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
        '''sampler of the TheAnnotator.

        Samples are created in a two step procedure. First, the length
        of a sample segment is chosen randomly from the empirical segment
        length distribution. Then, a random coordinate is chosen. If the
        sampled segment does not overlap with the workspace it is rejected.

        Before adding the segment to the sampled set, it is truncated to
        the workspace.

        If it overlaps with a previously sampled segment, the segments
        are merged. Thus, bases shared between two segments are not counted
        twice.

        Sampling continues, until exactly the same number of bases overlap between
        the ``P`` and the ``W`` as do between ``S`` and ``W``.

        Note that the length distribution of the ``P`` is different from ``S``.

        This method is quick if the workspace is large compared to the segments that
        need to be placed. If it is small, a large number of samples will be rejected.
        '''

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.nunsuccessful_rounds

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''

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
                start, end, overlap = temp_sampler.sample( length )
                unintersected_segments.trim( start, -true_remaining )

            # sample a position until we get a nonzero overlap
            start, end, overlap = sls.sample( length )
            sampled_segment = Segment( start, end )

            # If intersect > remaining we could well add too much (but we're not sure, since
            # there may be overlap with previously sampled segments).  Make sure we don't
            # overshoot the mark.  However, if we can't increase the amount of overlap in a
            # reasonable amount of time, don't bother
            if remaining < overlap and nunsuccessful_rounds < max_unsuccessful_rounds:
                # print "truncating sampled segment", remaining, overlap
                if workspace.overlapWithRange( pos, pos+1 ) > 0:
                    # Left end overlaps with workspace
                    sampled_segment.end = sampled_segment.start + remaining;
                else:
                    # Right end does?
                    sampled_segment.start = sampled_segment.end - remaining;
                overlap = workspace.overlapWithRange( sampled_segment.start, sampled_segment.end );

            if true_remaining > 0:
                # print "adding segment ", sampled_segment, true_remaining
                sampled_segments._add( sampled_segment )
                remaining -= overlap

        self.nunsuccessful_rounds = nunsuccessful_rounds
        return intersected_segments

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

    def __call__(self, SegmentList segments, SegmentList annotations, SegmentList workspace = None ):
        '''return number of nucleotides overlapping between segments and annotations.'''
        return annotations.overlapWithSegments( segments )

cdef class CounterNucleotideDensity(Counter):

    def __call__(self, SegmentList segments, SegmentList annotations, SegmentList workspace ):
        '''return number of nucleotides overlapping between segments and annotations.
        divided by the size of the workspace
        '''
        cdef long l
        l = len(workspace)
        if l == 0: return 0
        return float(annotations.overlapWithSegments( segments )) / l

cdef class CounterSegmentOverlap(Counter):

    def __call__(self, segments, annotations, workspace = None ):
        '''return number of segments overlapping with annotations.'''
        return annotations.intersectionWithSegments( segments )

############################################################
############################################################
############################################################
## Annotator results
############################################################
@cython.boundscheck(False)
cpdef getNPTwoSidedPValue( ar, val ):
    '''return pvalue for *val* within sorted array *ar*
    '''
    idx = numpy.searchsorted( ar, val )
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

    return dmax( min_pval, pval)

@cython.boundscheck(False)
def getNPTwoSidedPValueFast( numpy.ndarray[DTYPE_FLOAT_t, ndim=1] ar,
                           double val,
                           double mean ):
    '''return a two-sided pvalue.

    Fast if val is small or large.
    '''
    cdef long l, x
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

    return dmax( min_pval, pval )

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
    long nsamples
    double * samples
    int * sorted2sample
    int * sample2sorted
    double pvalue
    double qvalue
    int observed_idx

@cython.profile(True)
cdef inline int isSampleSignificantAtPvalue( EnrichmentStatistics * stats, 
                                             long sample_id, 
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
    cdef long l

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

cdef double getTwoSidedPValue( EnrichmentStatistics * stats, 
                                double val ):
    '''return pvalue for *val* within sorted array *ar*
    '''
    cdef long idx, l
    cdef double min_pval, pval, rval
    
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
    cdef long l
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

    stats.observed_idx = observed_idx
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

cdef EnrichmentStatistics * makeEnrichmentStatistics( observed, samples ):

    cdef EnrichmentStatistics * stats 
    cdef long offset, i, l

    l = len(samples)

    stats = <EnrichmentStatistics*>malloc( sizeof(EnrichmentStatistics) )
    stats.samples = <double*>calloc( l, sizeof(double))
    stats.sorted2sample = <int*>calloc( l, sizeof(int))
    stats.sample2sorted = <int*>calloc( l, sizeof(int))

    stats.observed = observed
    stats.nsamples = l 
    for i from 0 <= i < l: stats.samples[i] = samples[i]

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
    if stats.expected != 0:
        stats.fold = stats.observed / stats.expected
    else:
        stats.fold = 1.0

    stats.stddev = numpy.std(samples)

    offset = int(0.05 * l)
    stats.lower95 = stats.samples[stats.sorted2sample[ lmin( offset, l-1) ]]
    stats.upper95 = stats.samples[stats.sorted2sample[ lmax( l-offset, 0) ]]
    
    stats.pvalue = getTwoSidedPValue( stats, stats.observed )
    stats.qvalue = 1.0

    return stats

############################################################
############################################################
############################################################
## Annotator results
############################################################
cdef class AnnotatorResult(object):
    '''container for annotator results.'''

    format_observed = "%i"
    format_expected = "%6.4f"
    format_fold = "%6.4f"
    format_pvalue = "%6.4e"

    headers = ["track", "annotation",
               "observed", "expected",
               "CI95low", "CI95high",
               "stddev",
               "fold",
               "pvalue",
               "qvalue" ]

    cdef:
        EnrichmentStatistics * stats
        str track, annotation

    def __init__( self,
                  track,
                  annotation,
                  observed,
                  samples ):

        self.track = track
        self.annotation = annotation
        
        self.stats = makeEnrichmentStatistics( observed, samples )

    def __str__(self):

        if self.stats.nsamples < 10**6:
            format_pvalue = "%7.6f"
        else:
            format_pvalue = "%7.6e"

        return "\t".join( (self.track,
                           self.annotation,
                           self.format_observed % self.stats.observed,
                           self.format_expected % self.stats.expected,
                           self.format_expected % self.stats.lower95,
                           self.format_expected % self.stats.upper95,
                           self.format_expected % self.stats.stddev,
                           self.format_fold % self.stats.fold,
                           format_pvalue % self.stats.pvalue,
                           format_pvalue % self.stats.qvalue,
                           ) )

    def __dealloc__(self):
        free( self.stats.samples )
        free( self.stats.sorted2sample )
        free( self.stats.sample2sorted )
        free( self.stats )

    property track:
        def __get__(self): return self.track

    property annotation:
        def __get__(self): return self.annotation

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

    property qvalue:
        def __get__(self): return self.stats.qvalue

    property samples:
        def __get__(self): 
            cdef long x
            r = numpy.zeros( self.nsamples, dtype = numpy.float )
            for x from 0 <= x < self.nsamples:
                r[x] = self.samples[x]
            return r
    def isSampleSignificantAtPvalue( self, sample_id, double pvalue ):
        return isSampleSignificantAtPvalue( self.stats, sample_id, pvalue )

    def getSample( self, sample_id ):
        return self.stats.samples[sample_id]

############################################################
############################################################
############################################################
## compute counts
############################################################
def _append( counts, x,y,z): counts[y].append(z)
def _set( counts, x,y,z) : counts[x].__setitem__( y, z )
def _defdictfloat(): return collections.defaultdict(float)

def computeCounts( counter, aggregator, 
                   segments, annotations, workspace,
                   append = False ):
    '''collect counts from *counter* between all combinations of *segments* and *annotations*.

    *aggregator* determines how values are combined across isochores.

    If *append* is set, values for each track in *segments* are appended to a list for
    each track in *annotations*. This is useful for samples.

    Otherwise, a nested dictionary is returned.
    '''

    if append:
        counts = collections.defaultdict( list )
        f = _append
    else:
        counts = collections.defaultdict( _defdictfloat )
        f = _set

    isochores = workspace.keys()
    l = len(isochores)
    # collect counts per isochore
    for track in segments.tracks:
        segs = segments[track]
        for annotation in annotations.tracks:
            annos = annotations[annotation]
            vals = [ counter( segs[isochore], annos[isochore], workspace[isochore] ) for isochore in workspace.keys() ]
            # vals = numpy.array( [ counter( segs[isochore], annos[isochore], workspace[isochore] ) for isochore in workspace.keys() ],
            #                     dtype = numpy.float )
            # note: not truly iterator - generators currently not available in cython
            # vals = numpy.fromiter( \
            #     [ counter( segs[isochore], annos[isochore], workspace[isochore] ) for isochore in isochores ],
            #     dtype = numpy.float,
            #     count = l)
            f(counts, track, annotation, aggregator(vals) )

    return counts

############################################################
############################################################
############################################################
## FDR computation
############################################################
cdef double computeFalsePositiveRate( EnrichmentStatistics ** allstats,
                                      long nresults,
                                      long nsamples,
                                      double pvalue ):
    '''return the number of expected number of false positives
    for less than or equal to *pvalue*.
    '''

    cdef long y, nsample, nfp, total_nfp
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

cpdef computeFDR( annotator_results ):
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

    cdef long nresults
    cdef AnnotatorResult r
    cdef EnrichmentStatistics * r1
    cdef EnrichmentStatistics ** allstats
    cdef long nsample, nfp, R, x,y, total_nfp
    cdef long * nfps

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

    def __cinit__(self ):

        self.data = NULL
        self.fields = NULL
        self.index = 0

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.data = buffer
        self.update( buffer, nbytes )

    cdef present( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        self.update( buffer, nbytes )

    cdef copy( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        cdef int s
        s = sizeof(char) * nbytes
        self.data = <char*>malloc( s )
        memcpy( <char*>self.data, buffer, s )
        self.update( self.data, nbytes )

    cdef update( self, char * buffer, size_t nbytes ):
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

    def __getitem__( self, key ):

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

cdef class BedProxy( TupleProxy ):

   cdef track
   def __init__(self, track ):
       self.track = track

   property contig:
       def __get__(self):
           assert 0 < self.nfields
           return self.fields[0]

   property start:
       def __get__(self):
           assert 1 < self.nfields
           return atol(self.fields[1])

   property end:
       def __get__(self):
           assert 2 < self.nfields
           return atol(self.fields[2])

   property name:
       def __get__(self):
           assert 3 < self.nfields
           return self.fields[3]

   property track:
       def __get__(self):
           return self.track

class Track(object):
    '''bed track information.'''
    def __init__(self, line ):
        r= re.compile('([^ =]+) *= *("[^"]*"|[^ ]*)')

        self._d = {}
        for k, v in r.findall(line[:-1]):
            if v[:1]=='"':
                self._d[k] = v[1:-1]
            else:
                self._d[k] = v

        self._line = line[:-1]
    def __str__(self):
        return self._line

    def __getitem__(self, key): return self._d[key]
    def __setitem__(self, key,val): self._d[key] = val

class tsv_iterator:
    '''iterate over ``infile``.

    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile ):

        self.infile = infile

    def __iter__(self):
        return self

    def preparse( self, line ):
        return True

    def create( self):
        return TupleProxy()

    def next(self):

        cdef char * b, * cpy

        cdef TupleProxy r
        cdef size_t nbytes

        track = None

        while 1:

            line = self.infile.readline()
            if not line: break

            b = line
            # nbytes includes new line but not \0
            nbytes = len( line )

            # skip comments
            if (b[0] == '#'): continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n': continue

            if not self.preparse(line): continue

            # make sure that entry is complete
            if b[nbytes-1] != '\n':
                raise ValueError( "incomplete line at %s" % line )

            # chop off newline
            b[nbytes-1] = '\0'

            # create a copy
            cpy = <char*>malloc( nbytes )
            memcpy( cpy, b, nbytes )

            r = self.create()
            r.take( cpy, nbytes )
            return r

        raise StopIteration

class bed_iterator(tsv_iterator):
    '''iterate over ``infile``.

    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile ):
        tsv_iterator.__init__(self, infile )
        self.track = None

    def preparse(self, line ):
        if line.startswith("track"):
            self.track = Track( line )
            return False
        return True

    def create( self ):
        return BedProxy( self.track )

def _genie():
    return collections.defaultdict(SegmentList)

def readFromBed( filenames ):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *contig* and *track*.

    If no track is given, the *name* attribute is taken.
    '''
    cdef SegmentList l
    cdef BedProxy bed

    segment_lists = collections.defaultdict( _genie )

    for filename in filenames:
        infile = IOTools.openFile( filename, "r")
        for bed in bed_iterator( infile ):
            if bed.track:
                name = bed.track["name"]
            else:
                name = bed.name

            l = segment_lists[name][bed.contig]
            l.add( atol(bed.fields[1]), atol(bed.fields[2]) )

    return segment_lists

#####################################################################
#####################################################################
#####################################################################
## samples
#####################################################################
class IntervalCollection(object):
    '''a collection of intervals.

    Intervals (objects of type :class:`SegmentList`) are collected by track and isochore.
    '''

    def __init__(self, name ):
        self.intervals = collections.defaultdict( _genie )
        self.name = name

    def load( self, filenames ):
        '''load segments from filenames and pre-process them.'''
        self.intervals = readFromBed( filenames )

    def normalize( self ):
        '''normalize segment lists individually.

        Remove empty contigs.
        '''

        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                if len(segmentlist) == 0: 
                    del vv[contig]
                else:
                    segmentlist.normalize()

    def outputStats( self, outfile ):
        '''output segment statistics.'''

        outfile.write( "section\ttrack\tcontig\tnsegments\tlength\n" )

        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % len(segmentlist), 
                         "%i" % segmentlist.sum() ) ) + "\n" )

    def merge( self ):
        '''merge all tracks into a single segment list
        creating a new track called 'merged`
        '''
        merged = collections.defaultdict( SegmentList )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                merged[contig].extend( segmentlist )
        self.intervals["merged"] = merged
        self.normalize()

    def sum( self ):
        '''remove all intervals not overlapping with intervals in other.'''
        s = 0
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                s += segmentlist.sum()
        return s

    def prune( self, other ):
        '''remove all intervals not overlapping with intervals in other.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.intersect( other[contig] )

    def sort( self ):
        '''sort all intervals lists.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.sort()

    def check( self ):
        '''check all intervals lists.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                segmentlist.sort()

    def restrict( self, restrict ):
        '''remove all tracks except those in restrict.'''
        if restrict in (list, tuple, set):
            r = set(restrict)
        else:
            r = set( [restrict,] )

        for track in self.intervals.keys():
            if track not in r:
                del self.intervals[track]

    def toIsochores( self, isochores ):
        '''split contigs in each track into isochores.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.items():
                for other_track, other_vv in isochores.iteritems():
                    newlist = SegmentList( clone = segmentlist )
                    newlist.intersect( other_vv[contig] )
                    isochore = "%s.%s" % (contig, other_track)
                    vv[isochore] = newlist
                del vv[contig]

    def dump( self, outfile ):
        '''dump in bed format.'''

        for track, vv in self.intervals.iteritems():
            outfile.write("track name=%s\n" % track )
            for contig, segmentlist in vv.items():
                for start, end in segmentlist:
                    outfile.write( "%s\t%i\t%i\n" % (contig, start, end))

    @property
    def tracks(self): return self.intervals.keys()

    def __len__(self): return len(self.intervals)

    def keys(self): return self.intervals.keys()

    def add( self, track, contig, segmentlist ):
        self.intervals[track][contig] = segmentlist

    def __getitem__(self, key ):
        return self.intervals[key]

    def __contains__(self, key ):
        return key in self.intervals

    def iteritems(self):
        return self.intervals.iteritems()
    
    def outputOverlapStats( self, outfile, other ):

        outfile.write( "section\ttrack\tcontig\toverlap\tlength\tdensity\n" )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                length = segmentlist.sum()
                if length == 0: continue
                overlap = segmentlist.overlapWithSegments(other[contig])
                outfile.write( "\t".join( \
                        (self.name, track, contig, 
                         "%i" % overlap,
                         "%i" % length,
                         "%f" % (float(overlap) / length)) ) + "\n" )

#####################################################################
#####################################################################
#####################################################################
## samples
#####################################################################
cdef class Samples( object ):
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

    def __getitem__(self, track ):
        '''return all samples for track (as an :class:`IntervalCollection`)'''
        return self.samples[track]

    def load( self, track, sample_id, isochore ):
        raise ValueError( "loading called for uncached data" )

    def __delitem__(self, key ):
        del self.samples[key]

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

        ckey = <char*>calloc( sizeof(char) * 256, 1 )
        x = 0
        while not feof( findex ):
            fread( &keylen, sizeof( char), 1, findex )
            if feof(findex): break
            fread( ckey, sizeof( char), keylen, findex )
            fread( &pos, sizeof(off_t), 1, findex )
            # converts to persistent python object
            self.index[ckey] = pos
            
        fclose(findex)
        free(ckey)

        E.debug( "loaded index from %s: %i items" % (self.filename,len(self.index) ))

    def add( self, track, sample_id, isochore, segmentlist ):
        '''add a new *sample* for *track* and *isochore*, giving it *sample_id*.'''

        Samples.add( self, track, sample_id, isochore, segmentlist )

        cdef SegmentList seglist
        seglist = segmentlist
        cdef off_t pos
        cdef char * key
        cdef char keylen

        l = len(seglist)
        if l == 0: return

        pos = ftello( self.fcache )

        # cache structure is:
        # 1 * sizeof(unsigned char) - key length (max 255 chars)
        # 1 * sizeof(long) - number of segments (nsegments)
        # nsegments * sizeof(Segment) - the segment list
        tempkey = "%s-%i-%s" % (track,sample_id,isochore)
        assert len(tempkey) <= 255
        key = tempkey
        keylen = strlen( key ) + 1

        self.index[ key ] = pos
        fwrite( &seglist.nsegments, sizeof( long ), 1, self.fcache )
        toCompressedFile( <unsigned char *>seglist.segments,
                          sizeof( Segment) * seglist.nsegments,
                          self.fcache )

        # write to index
        fwrite( &keylen, sizeof(char), 1, self.findex )
        fwrite( key, sizeof(char), keylen, self.findex )
        fwrite( &pos, sizeof(off_t), 1, self.findex )

    def hasSample( self, track, sample_id, isochore ):
        '''return true if cache has sample.'''
        key = "%s-%i-%s" % (track,sample_id,isochore)
        return key in self.index

    def __dealloc__(self):
        fclose( self.fcache )
        fclose( self.findex )

    def __getitem__(self, track ):
        '''return all samples for track (as an :class:`IntervalCollection`)'''
        return self.samples[track]

    def load( self, track, sample_id, isochore ):
        '''load data into memory'''
        cdef off_t pos
        cdef SegmentList seglist
        cdef long nsegments

        tempkey = "%s-%i-%s" % (track,sample_id,isochore)
        pos = self.index[tempkey]
        fseeko( self.fcache, pos, SEEK_SET )
        fread( &nsegments, sizeof(long), 1, self.fcache )
        seglist = SegmentList( allocate = nsegments )
        seglist.nsegments = nsegments
        fromCompressedFile( <unsigned char*> seglist.segments,
                             sizeof( Segment) * seglist.nsegments,
                             self.fcache )

        Samples.add( self, track, sample_id, isochore, seglist )

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

#         for track, samples in self.samples.iteritems():
#             group = self.cache.createGroup( self.cache.root, track, track)
#             for sample_id, sample in samples.iteritems():
#                 subgroup = self.cache.createGroup( group, str(sample_id), str(sample_id) )
#                 for isochore, seglist in sample.iteritems():
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

