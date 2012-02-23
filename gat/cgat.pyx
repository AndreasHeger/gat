# cython: embedsignature=True
# cython: profile=False

import types, collections, re, os, math
import gat.IOTools as IOTools
import gat.Experiment as E
import gat.Stats

cimport cython

cdef extern from "string.h":
    ctypedef int size_t
    void *memcpy(void *dest, void *src, size_t n)
    void *memmove(void *dest, void *src, size_t n)
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

cdef extern from "math.h":
    double floor(double x)

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

#####################################################
## scipy import
## required for pvalues based on parametric distributions
try:
    import scipy
    import scipy.stats
    HAS_SCIPY = True
except:
    HAS_SCIPY = False

# type for positions
ctypedef unsigned int Position
# type for position difference - can be negative
ctypedef int PositionDifference

cdef struct Segment:
    Position start
    Position end

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
    cdef Segment * segments
    cdef size_t nsegments
    cdef size_t allocated
    cdef int is_normalized
    cdef int chunk_size

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

    cpdef trim_ends( self, Position pos, Position size, int forward = 1 ):
        '''trim segment list by removing *size* nucleotides from
        the segment that includes *pos*.

        *forward* gives the direction.
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

    cpdef merge( self, PositionDifference distance = 0 ):
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

    cpdef Position overlapWithRange( self, Position start, Position end ):
        '''return the size of intersection between
           segment list and Segment other'''

        cdef Segment s
        s = Segment( start, end )
        return self.overlap( s )

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

    cpdef Position intersectionWithSegments( self, SegmentList other ):
        '''return the number of segments overlapping with *other*.'''

        assert self.is_normalized, "intersection from non-normalized list"
        assert other.is_normalized, "intersection with non-normalized list"

        # avoid self-self comparison
        if other.segments == self.segments: return self.sum()

        cdef int this_idx, other_idx, working_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Segment this_segment = Segment(0,0)
        cdef Segment other_segment = Segment(0,0)

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

cdef class SegmentListSamplerSlow:

    cdef SegmentList segment_list
    cdef numpy.ndarray  cdf
    cdef Position total_size

    def __init__(self, SegmentList segment_list ):

        assert len(segment_list) > 0, "sampling from empty segment list"

        self.segment_list = segment_list
        self.cdf = numpy.cumsum( [x[1] - x[0] for x in self.segment_list ] )
        self.total_size = self.cdf[-1]

    cpdef sample( self, Position length ):
        '''return a new position within segment list.'''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef Position r, offset, pos
        cdef PositionDifference overlap
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
    cdef Position * cdf
    cdef Position total_size
    cdef Position nsegments

    def __init__(self, SegmentList segment_list ):
        cdef Position i, totsize
        assert len(segment_list) > 0, "sampling from empty segment list"

        self.segment_list = segment_list
        self.nsegments = len(segment_list)
        self.cdf = <Position*>malloc( sizeof(Position) * self.nsegments )
        self.total_size = 0
        for i from 0 <= i < len(segment_list):
            self.total_size += segment_length( segment_list.segments[i] )
            self.cdf[i] = self.total_size

    cpdef sample( self, Position length ):
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

        rseg = numpy.random.randint( 0, length )

        start, end = pos - rseg, pos - rseg + length

        # print "sample %i-%i, offset=%i, pos=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
        #     (start, end, offset, pos, rpos, rseg,
        #      segment_index, self.segment_list.nsegments,
        #      self.segment_list.segments[segment_index] )
        overlap = 1

        assert overlap > 0, "sample %i-%i does not overlap with workspace: offset=%i, rpos=%i, rseg=%i, index=%i/%i, segment=%s" %\
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
        assert segment_list.is_normalized

        self.segment_list = segment_list
        self.nsegments = len(segment_list)
        self.cdf = <Position*>malloc( sizeof(Position) * self.nsegments )
        self.total_size = 0
        for i from 0 <= i < len(segment_list):
            self.total_size += segment_length( segment_list.segments[i] )
            # -1, as searchsorted inserts on left
            self.cdf[i] = self.total_size -1

    cpdef sample( self, Position sample_length ):
        '''return a sample.'''

        # note: could be made quicker by
        # - creating a sample random integers once, see numpy.random_integers?
        # - avoiding the binary search?
        cdef Position start, end
        cdef PositionDifference overlap
        cdef size_t segment_index
        cdef Segment chosen_segment

        cdef Position random_pos_in_workspace, random_pos_in_segment

        # sample a position 
        # r = rand() / (RAND_MAX / total_size + 1)
        random_pos_in_workspace = numpy.random.randint( 0, self.total_size )
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
        random_pos_in_segment = numpy.random.randint( lmax(0, chosen_segment.start - sample_length + 1),
                                                      chosen_segment.end )

        # adjust position by sampled amount
        start = random_pos_in_segment
        end = start + sample_length

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
                 Position bucket_size ):

        assert len(histogram) > 0, "sampling from empty histogram"

        self.cdf = numpy.cumsum( histogram, dtype = numpy.float )
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
            free( self.cdf )
            self.cdf = NULL

cdef class Sampler:
    pass

cdef class SamplerAnnotator(Sampler):

    cdef Position bucket_size
    cdef int nbuckets
    cdef int nunsuccessful_rounds

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

        assert segments.is_normalized, "segment list is not normalized"
        assert workspace.is_normalized, "workspace is not normalized"

        unintersected_segments = SegmentList()
        intersected_segments = SegmentList()

        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )

        # or: intersect?
        working_segments.filter( workspace )

        if len(working_segments) == 0:
            return intersected_segments

        # get nucleotides that need to be sampled
        # only count overlap within workspace
        tmp_segments = SegmentList( clone = working_segments )
        tmp_segments.intersect( workspace )
        ltotal = tmp_segments.sum()        

        # create space for sampled segments, add additional 10%
        # safety margin to avoid realloc calls
        sampled_segments = SegmentList( allocate = int(1.1 * len(segments)) )

        # build length histogram
        histogram = working_segments.getLengthDistribution( self.bucket_size,
                                                            self.nbuckets )

        hs = HistogramSampler( histogram, self.bucket_size )

        # get segment sampler
        sls = SegmentListSampler( workspace )

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
                #print "recomputing segments", nunsuccessful_rounds
                #print " sampled segments", sampled_segments.asList()
                unintersected_segments.extend( sampled_segments )
                unintersected_segments.merge()
                #print "extended"
                sampled_segments.clear()
                # compute overlap with workspace
                intersected_segments = SegmentList( clone = unintersected_segments )
                intersected_segments.intersect( workspace )

                #print " unintersected_segments", unintersected_segments.sum(), len(unintersected_segments)
                #print unintersected_segments.asList()
                #print " intersected_segments", intersected_segments.sum(), len(intersected_segments)
                #print intersected_segments.asList()
                assert intersected_segments.sum() <= unintersected_segments.sum()
                #if not intersected_segments.sum() <= ltotal:
                #   print "\n".join( [str(x) for x in intersected_segments ] ) + "\n"
                #assert intersected_segments.sum() <= ltotal, \
                #    "length of segments exceeds ltotal: %i > %i" % (intersected_segments.sum(), ltotal)
                remaining = ltotal - <PositionDifference>intersected_segments.sum()
                if true_remaining == remaining:
                    nunsuccessful_rounds += 1
                else:
                    true_remaining = remaining

            # deal with overshoot
            if true_remaining < 0:
                #print "removing overshoot", true_remaining
                temp_sampler = SegmentListSampler( unintersected_segments )
                # sample a position from current list of sampled segments
                start, end, overlap = temp_sampler.sample( 1 )
                #print "unintersected_segments before trim", unintersected_segments.sum(), unintersected_segments.asList()
                # causes a potential bias for overlapping segments- in case of overlapping segments, the
                # segment chosen is always the first
                unintersected_segments.trim_ends( start, 
                                                  -true_remaining,
                                                  numpy.random.randint( 0,2) )
                # trimming might remove residues outside the workspace
                # hence - recompute true_remaining by going back to loop.
                #print "unintersected_segments after  trim", unintersected_segments.sum(), unintersected_segments.asList()
                true_remaining = 1
                continue

            # sample a position until we get a nonzero overlap
            start, end, overlap = sls.sample( length )
            #print start, end, overlap, remaining, true_remaining
            sampled_segment = Segment( start, end )

            # If intersect > remaining we could well add too much (but we're not sure, since
            # there may be overlap with previously sampled segments).  Make sure we don't
            # overshoot the mark.  However, if we can't increase the amount of overlap in a
            # reasonable amount of time, don't bother
            # 
            #     # note that this part causes an edge bias as segments are preferentially added
            #     # from within a workspace, thus depleting counts at positions were the workspace 
            #     # is discontinuous.
            # 
            # if remaining < overlap and nunsuccessful_rounds < max_unsuccessful_rounds:
            #     #print "truncating sampled segment", remaining, overlap, str(sampled_segment)
            #     if workspace.overlapWithRange( sampled_segment.start, sampled_segment.start+1 ) > 0:
            #         # Left end overlaps with workspace
            #         sampled_segment.end = sampled_segment.start + remaining;
            #     elif workspace.overlapWithRange( sampled_segment.end-1, sampled_segment.end ) > 0:
            #         # Right end does?
            #         sampled_segment.start = sampled_segment.end - remaining;

            #     overlap = workspace.overlapWithRange( sampled_segment.start, sampled_segment.end );
            #     #print "truncated to ", str(sampled_segment)

            if true_remaining > 0:
                # print "adding segment ", sampled_segment, true_remaining
                sampled_segments._add( sampled_segment )
                remaining -= overlap

        self.nunsuccessful_rounds = nunsuccessful_rounds
        
        intersected_segments = SegmentList( clone = unintersected_segments )
        #print "final: unintersected", intersected_segments.sum(), intersected_segments.asList()
        intersected_segments.merge()
        #print "final: merged", intersected_segments.sum(), intersected_segments.asList()
        intersected_segments.filter( workspace )
        #print "final: filtered", intersected_segments.sum(), intersected_segments.asList()
        #print "final:", intersected_segments.sum(), intersected_segments.asList()
        
        assert intersected_segments.sum() > 0
        return intersected_segments

########################################################
########################################################
########################################################
cdef class SamplerSegments(Sampler):

    cdef Position bucket_size
    cdef Position nbuckets

    def __init__( self, bucket_size = 1, nbuckets = 100000 ):
        '''sample *n* segments from length distribution.

        The segment length distribution is derived from
        the argument *segments* supplied to the :meth:`sample`
        method.

        *n* is given by the number of segments in the observed
        data.

        *segments* does not need to be normalized. In fact,
        supplying overlapping segments is likely to be the 
        most realistic use case.
        '''

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''

        cdef Position length
        cdef SegmentListSampler sls
        cdef SegmentList sample

        assert workspace.is_normalized, "workspace is not normalized"

        sample = SegmentList( allocate = len(segments) )
    
        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )
        working_segments.filter( workspace )

        if len(working_segments) == 0:
            return sample

        # build length histogram
        histogram = working_segments.getLengthDistribution( self.bucket_size,
                                                            self.nbuckets )
        hs = HistogramSampler( histogram, self.bucket_size )

        # create segment sampler
        sls = SegmentListSampler( workspace )

        for x in xrange( len(segments) ):

            # Sample a segment length from the histogram
            length = hs.sample()
            assert length > 0

            # sample a position until we get a nonzero overlap
            start, end, overlap = sls.sample( length )

            sample._add( Segment( start, end ) )

        return sample

########################################################
########################################################
########################################################
cdef class SamplerBruteForce(Sampler):

    cdef Position bucket_size
    cdef Position nbuckets
    cdef int ntries_inner
    cdef int ntries_outer

    def __init__( self, 
                  bucket_size = 1, 
                  nbuckets = 100000,
                  ntries_inner = 100,
                  ntries_outer = 10 ):
        '''sample segments from length distribution until
        workspace is covered by the same number of nucleotides
        as input.

        The segment length distribution is derived from
        the argument *segments* supplied to the :meth:`sample`
        method.

        Sample by brute force - add and remove segments until
        the correct number of nucleotides has been sampled.
        '''

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.ntries_inner = ntries_inner
        self.ntries_outer = ntries_outer

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''

        cdef Position start, end, length
        cdef SegmentListSampler sls
        cdef SegmentList sample
    
        assert workspace.is_normalized, "workspace is not normalized"

        sample = SegmentList( allocate = len(segments) )

        # collect all segments in workspace
        working_segments = SegmentList( clone = segments )
        working_segments.filter( workspace )

        if len(working_segments) == 0:
            return sample

        # build length histogram
        histogram = working_segments.getLengthDistribution( self.bucket_size,
                                                            self.nbuckets )
        hs = HistogramSampler( histogram, self.bucket_size )

        # create segment sampler
        sls = SegmentListSampler( workspace )

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
                start, end, overlap = sls.sample( length )

                # print str(sample), start, end, overlap, remaining

                if overlap > remaining:
                    ntries_inner -= 1
                    continue

                if sample.overlapWithRange( start, end ): 
                    ntries_inner -= 1
                    continue

                sample._add( Segment( start, end ) )
                # required in order to sort samples
                sample.normalize()

                # print "adding", start, end, str(sample)

                ntries_inner = self.ntries_inner

                remaining -= overlap

            if ntries_inner > 0:
                break

            ntries_outer -= 1
            
        if ntries_outer == 0:
            raise ValueError( "sampling did not converge: %s" % str(sample))

        return sample

########################################################
########################################################
########################################################
cdef class SamplerUniform(Sampler):

    cdef Position bucket_size
    cdef Position nbuckets
    cdef Position increment
    cdef Position start_at
    cdef Position current_workspace
    cdef Position current_position
    cdef int current_orientation

    def __init__( self, 
                  increment,
                  bucket_size = 1,
                  nbuckets = 100000 ):
        '''
        For debugging purposes.

        Every *increment* residues within the worspace
        a new segment is generated. The segment length is
        taken from the segment size distribution.

        Segments are alternately forward/backward looking.
        '''

        self.bucket_size = bucket_size
        self.nbuckets = nbuckets
        self.increment = increment
        self.current_orientation = 0
        self.current_workspace = 0
        self.current_position = 0

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''
        assert workspace.is_normalized, "workspace is not normalized"

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
        histogram = working_segments.getLengthDistribution( self.bucket_size,
                                                            self.nbuckets )

        hs = HistogramSampler( histogram, self.bucket_size )

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

########################################################
########################################################
########################################################
cdef class SamplerShift(Sampler):

    cdef double radius

    def __init__( self, 
                  radius = 2 ):
        '''
        Sample segments by shifting them by a random amount within *radius*
        in a random direction.

        *radius* determines the size of the region that is accessible 
        for randomly shifting a segment. It is expressed as a fraction 
        of the size of a segment.
        
        '''

        self.radius = radius

    cpdef SegmentList sample( self,
                              SegmentList segments,
                              SegmentList workspace ):
        '''return a sampled list of segments.'''
        assert workspace.is_normalized, "workspace is not normalized"

        cdef Segment segment
        cdef SegmentList sample, working_segments
        cdef Position length, direction, extended_length
        cdef Position x 
        cdef PositionDifference shift, start, end
        cdef double half_radius = self.radius / 2
        cdef Segment * _working_segments

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
            shift_area = <Position>floor(length * half_radius)
            shift = numpy.random.randint( -shift_area, shift_area  )
            start = lmax(0, segment.start + shift )
            end = lmax( 0, segment.end + shift )
            sample._add( Segment(start, end) )

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

    def __call__(self, SegmentList segments, SegmentList annotations, SegmentList workspace = None ):
        '''return number of nucleotides overlapping between segments and annotations.'''
        return annotations.overlapWithSegments( segments )

cdef class CounterNucleotideDensity(Counter):

    def __call__(self, SegmentList segments, SegmentList annotations, SegmentList workspace ):
        '''return number of nucleotides overlapping between segments and annotations.
        divided by the size of the workspace
        '''
        cdef Position l
        l = len(workspace)
        if l == 0: return 0
        return float(annotations.overlapWithSegments( segments )) / l

cdef class CounterSegmentOverlap(Counter):

    def __call__(self, segments, annotations, workspace = None ):
        '''return number of segments overlapping with annotations.'''
        return segments.intersectionWithSegments( annotations )

cdef class CounterAnnotationOverlap(Counter):

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

cdef EnrichmentStatistics * makeEnrichmentStatistics( observed, samples ):

    cdef EnrichmentStatistics * stats 
    cdef Position offset, i, l

    l = len(samples)
    if l < 1:
        raise ValueError( "not enough samples - no stats in makeEnrichmentStatistics" )

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
    
    if offset > 0: 
        stats.lower95 = stats.samples[stats.sorted2sample[ lmin( offset, l-1) ]]
        stats.upper95 = stats.samples[stats.sorted2sample[ lmax( l-offset, 0) ]]
    else:
        stats.lower95 = stats.samples[stats.sorted2sample[ 0 ]]
        stats.upper95 = stats.samples[stats.sorted2sample[ l-1 ]]
    
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
        EnrichmentStatistics * stats
        str track, annotation
        Position track_nsegments
        Position track_size
        Position annotation_nsegments
        Position annotation_size
        Position overlap_nsegments
        Position overlap_size
        Position workspace_size

    def __init__( self,
                  track,
                  annotation,
                  observed,
                  samples,
                  track_segments,
                  annotation_segments,
                  workspace ):

        self.track = track
        self.annotation = annotation
        self.stats = makeEnrichmentStatistics( observed, samples )

        self.track_nsegments = track_segments.counts()
        self.track_size = track_segments.sum()

        self.annotation_nsegments = annotation_segments.counts()
        self.annotation_size = annotation_segments.sum()

        overlap = track_segments.clone()
        overlap.intersect( annotation_segments )

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
                           self.format_density % ( float(self.track_size) / self.workspace_size),
                           self.format_counts % self.annotation_nsegments,
                           self.format_counts % self.annotation_size,
                           self.format_density % ( float(self.annotation_size) / self.workspace_size),
                           self.format_counts % self.overlap_nsegments,
                           self.format_counts % self.overlap_size,
                           self.format_density % ( float(self.overlap_size) / self.workspace_size),
                           self.format_fold % ( 100.0 * float( self.overlap_nsegments) / self.track_nsegments ),
                           self.format_fold % ( 100.0 * float( self.overlap_size) / self.track_size ),
                           self.format_fold % ( 100.0 * float( self.overlap_nsegments) / self.annotation_nsegments ),
                           self.format_fold % ( 100.0 * float( self.overlap_size) / self.annotation_size ),
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

    def isSampleSignificantAtPvalue( self, sample_id, double pvalue ):
        return isSampleSignificantAtPvalue( self.stats, sample_id, pvalue )

    def getSample( self, sample_id ):
        return self.stats.samples[sample_id]

    def getEmpiricalPValue( self, value ):
        return getTwoSidedPValue( self.stats, value )

def getNormedPValue( value, r ):
    '''return pvalue assuming that samples are normal distributed.'''
    absval = abs(value - r.expected)
    if HAS_SCIPY:
        pvalue = 1.0 - scipy.stats.norm.cdf( absval, 0, r.stddev )
    else:
        raise ImportError( "scipy required" )
    return pvalue

def getEmpiricalPValue( value, r ):
    return r.getEmpiricalPValue( value )

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

def updateQValues( annotator_results, method = "storey", **kwargs ):
    '''update qvalues in annotator results

    storey
        qvalue from the method by Storey et al.
        print track, other_track, contig
    '''

    pvalues = [ r.pvalue for r in annotator_results ]

    if method == "storey":
        try:
            fdr = gat.Stats.computeQValues( pvalues, 
                                            vlambda = kwargs.get( "vlambda", numpy.arange( 0,0.95,0.05) ),
                                            pi0_method = kwargs.get( "pi0_method", "smoother") )
        except ValueError, msg:
            E.warn( "qvalue computation failed: %s" % msg )
            return

        for r, qvalue in zip(annotator_results, fdr.qvalues):
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
            temp_segments.filter( workspace )
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

class ConditionalWorkspaceCentered( UnconditionalWorkspace ):
    '''a workspace centered around segments/annotations.'''

    is_conditional = True

    def __init__( self, extension = None, expansion = None ):
        self.extension = extension
        self.expansion = expansion
        if self.extension == self.expansion == None:
            raise ValueError( "need to specify either expansion or extension" )

    def __call__( self, segments, annotations, workspace ):
        
        # build workspace from annotations
        temp_workspace = self.getCenter( segments, annotations ).clone()
        if self.extension != None:
            temp_workspace.extend( self.extension )
        else:
            temp_workspace.expand( self.expansion )
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

def computeCounts( counter, 
                   aggregator, 
                   segments, 
                   annotations, 
                   workspace,
                   workspace_generator,
                   append = False):
    '''collect counts from *counter* between all combinations of *segments* and *annotations*.

    *aggregator* determines how values are combined across isochores.

    If *append* is set, values for each track in *segments* are appended to a list for
    each track in *annotations*. This is useful for samples.  Otherwise, a nested dictionary is returned.

    '''
        
    if append:
        counts = collections.defaultdict( list )
        f = _append
    else:
        counts = collections.defaultdict( _defdictfloat )
        f = _set

    isochores = workspace.keys()

    # collect counts per isochore and aggregate these
    for track in segments.tracks:
        segs = segments[track]
            
        for annotation in annotations.tracks:
            annos = annotations[annotation]

            temp_segs, temp_annos, temp_workspace = workspace_generator( segs, annos, workspace )
            vals = [ counter( segs[isochore], annos[isochore], workspace[isochore] ) for isochore in isochores ]
            f(counts, track, annotation, aggregator(vals) )

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

    cdef copy( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        if self.data != NULL:
            free(self.data)

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

        if self.nfields < self.min_fields:
            raise IOError("not enough fields, expected >%i, got %i" % \
                              (self.min_fields, self.nfields))

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

cdef class BedProxy( TupleProxy ):

   cdef track
   def __init__(self, track ):
       self.track = track
       self.min_fields = 3

   property contig:
       def __get__(self):
           return self.fields[0]

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
    def __init__(self, line ):
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

        cdef char * b

        cdef TupleProxy r
        cdef size_t nbytes

        track = None

        while 1:

            line = self.infile.readline()
            if not line: break

            # string conversion - b references internal string in python object
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
            r = self.create()
            r.copy( b, nbytes )
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
    return IntervalCollection(SegmentList)

def readFromBed( filenames, allow_multiple = False ):
    '''read Segment Lists from one or more bed files.

    Segment lists are grouped by *contig* and *track*.

    If no *track* is given, the *name* attribute is taken
    instead. If that is not present (BED3 format), the
    *track* is called ``default``.

    If a *track* appears in multiple files, an error is raised
    unless *allow_multiple* is ``True``.
    
    '''
    cdef SegmentList l
    cdef BedProxy bed
    cdef Position lineno

    segment_lists = collections.defaultdict( IntervalDictionary )

    tracks = {}

    if type(filenames) == str:
        filenames = (filenames,)
    for filename in filenames:
        infile = IOTools.openFile( filename, "r")
        default_name = os.path.basename( filename )
        lineno = 0
        try:
            for bed in bed_iterator( infile ):
                if bed.track:
                    try:
                        name = bed.track["name"]
                    except KeyError:
                        raise KeyError( "track without field 'name' in file '%s'" % filename)
                elif bed.name: 
                    name = bed.name
                else:
                    name = default_name

                if name in tracks:
                    if tracks[name] != filename: 
                        if allow_multiple:
                            E.warn( "track '%s' in multiple filenames: %s and %s" % (name, tracks[name], filename))
                        else:
                            raise ValueError( "track '%s' in multiple filenames: %s and %s" % (name, tracks[name], filename) )
                        tracks[name] = filename
                else:
                    tracks[name] = filename
                 
                l = segment_lists[name][bed.contig]
                l.add( atol(bed.fields[1]), atol(bed.fields[2]) )
                lineno += 1

        except IOError, msg:
            raise IOError("malformatted entry in line %s:%i, msg=%s" % (filename,
                                                                        lineno,
                                                                        msg))

    return segment_lists

#####################################################################
#####################################################################
#####################################################################
## 
#####################################################################
class IntervalDictionary( object ):
    '''a collection of intervals.
    
    Intervals (objects of type :class:`SegmentList`) are organized
    in hierarchical dictionary by isochore.
    '''

    def __init__(self ):
        self.intervals = collections.defaultdict( SegmentList )

    def sum( self ):
        '''return sum of all segment lists.'''
        s = 0
        for contig, segmentlist in self.intervals.iteritems():
            s += segmentlist.sum()
        return s

    def counts( self ):
        '''return number of all segments in all segments lists.'''
        s = 0
        for contig, segmentlist in self.intervals.iteritems():
            s += len(segmentlist)
        return s

    def intersect( self, other ):
        '''intersect with intervals in other.'''
        for contig, segmentlist in self.intervals.items():
            if contig in other:
                segmentlist.intersect( other[contig] )
            else:
                del self.intervals[contig]

    def extend( self, extension ):
        '''extend each interval by a certain amount.'''
        for contig, segmentlist in self.intervals.items():
            segmentlist.extend_segments( extension )

    def expand( self, expansion ):
        '''expand each interval by a certain amount.'''
        for contig, segmentlist in self.intervals.items():
            segmentlist.expand_segments( expansion )

    def normalize( self ):
        '''normalize all segment lists.'''
        for contig, segmentlist in self.intervals.items():
            segmentlist.normalize()

    def __len__(self): return len(self.intervals)

    def __delitem__(self,key): del self.intervals[key]

    def keys(self): return self.intervals.keys()

    def add( self, contig, segmentlist ):
        self.intervals[contig] = segmentlist

    def items(self ):
        return self.intervals.items()

    def __getitem__(self, key ):
        return self.intervals[key]

    def __setitem__(self, key, val ):
        self.intervals[key] = val

    def __contains__(self, key ):
        return key in self.intervals

    def iteritems(self):
        return self.intervals.iteritems()

    def clone( self ):
        '''return a copy of the data.'''
        r = IntervalDictionary()
        for contig, segmentlist in self.intervals.iteritems():
            r[contig] = SegmentList( clone = segmentlist )
        return r

    def filter( self, other ):
        '''remove all intervals not overlapping with intervals in other.'''
        for contig, segmentlist in self.intervals.items():
            if contig in other:
                segmentlist.filter( other[contig] )
            else:
                del self.intervals[contig]

#####################################################################
#####################################################################
#####################################################################
## 
#####################################################################
class IntervalCollection(object):
    '''a collection of intervals.

    Intervals (objects of type :class:`SegmentList`) are organized
    in hierarchical dictionary by track and isochore.
    '''

    def __init__(self, name ):
        self.intervals = collections.defaultdict( IntervalDictionary )
        self.name = name

    def load( self, filenames, split_tracks = False ):
        '''load segments from filenames and pre-process them.'''
        self.intervals = readFromBed( filenames, split_tracks )
        self.normalize()

    def save( self, outfile, prefix = "", **kwargs ):
        '''save in bed format to *outfile*.
        
        Each interval set will be saved as a different track with optional
        prefix *prefix*.

        Additional *kwargs* will be added to the tracks line as key=value pairs.
        '''

        for track, vv in self.intervals.iteritems():
            outfile.write("track name=%s%s %s\n" % \
                              (prefix, track,
                               " ".join( ["%s=%s" % (x,y) for x,y in kwargs.iteritems() ] ) ) )
            for contig, segmentlist in vv.items():
                for start, end in segmentlist:
                    outfile.write( "%s\t%i\t%i\n" % (contig, start, end))

    def normalize( self ):
        '''normalize segment lists individually.

        Remove empty contigs.
        '''

        for track, vv in self.intervals.iteritems():
            for contig in vv.keys():
                segmentlist = vv[contig]
                if len(segmentlist) == 0: 
                    del vv[contig]
                else:
                    segmentlist.normalize()

    def outputStats( self, outfile ):
        '''output segment statistics.'''

        outfile.write( "section\ttrack\tcontig\tnsegments\tlength\n" )

        cdef Position total_length = 0
        cdef Position total_segments = 0
        cdef Position length, segments
        for track, vv in self.intervals.iteritems():
            total_length, total_segments = 0, 0
            for contig, segmentlist in vv.iteritems():
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
                 
                
    def merge( self, delete = False ):
        '''merge all tracks into a single segment list
        creating a new track called 'merged'.

        Overlapping intervals will be merged.

        If *delete* == True, all tracks but the merged track
        will be deleted.
        '''
        merged = IntervalDictionary()
        
        for track in self.intervals.keys():
            vv = self.intervals[track]
            for contig, segmentlist in vv.iteritems():
                merged[contig].extend( segmentlist )
            if delete:
                del self.intervals[track]

        self.intervals["merged"] = merged
        self.normalize()

    def collapse( self ):
        '''collapse all tracks into a single segment list
        creating a new track called 'collapsed'.

        The collapsed track contains the intersection of
        all segment lists.
        '''
        result = IntervalDictionary()
        
        # get list of contigs present in all data sets
        contigs = collections.defaultdict( int )
        for track, vv in self.intervals.iteritems():
            for contig in vv.keys():
                contigs[contig] += 1

        ntracks = len( self.intervals )
        shared_contigs = set( [ x for x,y in contigs.iteritems() if y == ntracks ] )
        
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                if contig not in shared_contigs: continue
                if contig not in result:
                    result[contig] = SegmentList( clone = segmentlist )
                else:
                    result[contig].intersect( segmentlist )
                    
        self.intervals["collapsed"] = result

    def sum( self ):
        '''return sum of all segment lists.'''
        s = 0
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                s += segmentlist.sum()
        return s

    def counts( self ):
        '''return number of all segments in all segments lists.'''
        s = 0
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.iteritems():
                s += len(segmentlist)
        return s

    def countsPerTrack( self ):
        '''return number of all segments in all segments lists.'''
        counts = {}
        for track, vv in self.intervals.iteritems():
            s = 0
            for contig, segmentlist in vv.iteritems():
                s += len(segmentlist)
            counts[track] = s 
        return counts
        
    def intersect( self, other ):
        '''intersect with intervals in other.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.items():
                if contig in other:
                    segmentlist.intersect( other[contig] )
                else:
                    del vv[contig]

    def filter( self, other ):
        '''remove all intervals not overlapping with intervals in other.'''
        for track, vv in self.intervals.iteritems():
            for contig, segmentlist in vv.items():
                if contig in other:
                    segmentlist.filter( other[contig] )
                else:
                    del vv[contig]

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

    @property
    def tracks(self): return self.intervals.keys()

    def __len__(self): return len(self.intervals)

    def __delitem__(self,key): del self.intervals[key]

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
            s.load( filename )
            self.samples[track] = s
            
    def load( self, track, sample_id, isochore ):
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
        # 1 * sizeof(Position) - number of segments (nsegments)
        # nsegments * sizeof(Segment) - the segment list
        tempkey = self.toKey( track, sample_id, isochore )
        assert len(tempkey) <= 255
        key = tempkey
        keylen = strlen( key ) + 1

        self.index[ key ] = pos

        fwrite( &seglist.nsegments, sizeof( Position ), 1, self.fcache )
        toCompressedFile( <unsigned char *>seglist.segments,
                          sizeof( Segment) * seglist.nsegments,
                          self.fcache )

        # write to index
        fwrite( &keylen, sizeof(char), 1, self.findex )
        fwrite( key, sizeof(char), keylen, self.findex )
        fwrite( &pos, sizeof(off_t), 1, self.findex )

    def toKey( self, track, sample_id, isochore ):
        return "%s-%s-%s" % (track,sample_id,isochore) 

    def hasSample( self, track, sample_id, isochore ):
        '''return true if cache has sample.'''
        return self.toKey( track, sample_id, isochore ) in self.index

    def __dealloc__(self):
        fclose( self.fcache )
        fclose( self.findex )

    def load( self, track, sample_id, isochore ):
        '''load data into memory'''
        cdef off_t pos
        cdef SegmentList seglist
        cdef Position nsegments

        tempkey = self.toKey( track, sample_id, isochore )
        pos = self.index[tempkey]
        fseeko( self.fcache, pos, SEEK_SET )
        fread( &nsegments, sizeof(Position), 1, self.fcache )
        seglist = SegmentList( allocate = nsegments )
        seglist.nsegments = nsegments
        fromCompressedFile( <unsigned char*> seglist.segments,
                             sizeof( Segment) * seglist.nsegments,
                             self.fcache )

        Samples.add( self, track, sample_id, isochore, seglist )

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

