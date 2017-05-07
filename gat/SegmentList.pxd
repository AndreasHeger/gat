from libc.stdio cimport FILE
from posix.types cimport off_t

from CoordinateList cimport CoordinateList

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

    int toCompressedFile(unsigned char *, size_t, FILE *)
    int fromCompressedFile(unsigned char *, size_t, FILE *)

cdef bytes force_bytes(object s, encoding=*)
cdef force_str(object s, encoding=*)

#####################################################
#####################################################
## type definitions
#####################################################
# type for positions
ctypedef unsigned int Position
# type for position difference - can be negative
ctypedef int PositionDifference

# A segment (without strand)
cdef struct Segment:
    Position start
    Position end

#####################################################
#####################################################
## definition of segmentlist
#####################################################
cdef class SegmentList(CoordinateList):

    cdef Segment * segments
    cdef size_t nsegments
    cdef size_t allocated
    cdef int flag

    # C and Python methods
    cpdef share(self, key)
    cpdef sort(self)
    cpdef SegmentList extend(self, SegmentList other)
    cpdef add(self, Position start, Position end)
    cpdef trim_ends(self, Position pos, Position size, int forward )
    cpdef trim(self, Position pos, Position size)
    cpdef normalize(self)
    cpdef merge(self, PositionDifference distance )
    cpdef check(self)
    cpdef Position getRandomPosition(self)
    cpdef Position overlapWithRange(
        self, Position start, Position end)
    cpdef SegmentList getOverlappingSegmentsWithRange(
        self, Position start, Position end)
    cpdef Position overlapWithSegments(
        self, SegmentList other)
    cpdef Position intersectionWithSegments(
        self, SegmentList other, mode = *)
    cpdef SegmentList getFilledSegmentsFromStart(
        self, Position start, PositionDifference remainder)
    cpdef SegmentList getFilledSegmentsFromEnd(
        self, Position end, PositionDifference remainder)
    cpdef void filter(self, SegmentList other)
    cpdef void intersect(self, SegmentList other)
    cpdef void subtract(self, SegmentList other)
    cpdef extend_segments(self, int extension)
    cpdef expand_segments(self, double expansion)
    cpdef Position sum(self)
    cpdef Position max(self)
    cpdef Position min(self)
    cpdef clear(self)
    cpdef summarize(self)
    cpdef SegmentList clone(self)

    # C only methods
    cdef _add(self, Segment segment)
    cdef _resize(self, int nsegments)
    cdef insert(self, int idx, Segment seg)
    cdef int _getInsertionPoint(self, Segment other)
    cdef Position overlap(self, Segment other)
    cdef SegmentList getOverlappingSegments(self, Segment other)
    cdef void truncate(self, Segment other)

    # Functions overloaded from CoordinateList
    cdef off_t toMMAP(self, void *, int, off_t)
    cdef void fromMMAP(self)

