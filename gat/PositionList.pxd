from posix.types cimport off_t
from SegmentList cimport Position, Segment, SegmentList
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


cdef class PositionList(CoordinateList):

    cdef Position * positions
    cdef size_t npositions
    cdef size_t allocated
    cdef bint is_sorted	
    cdef bint is_normalized

    cpdef sort(self)
    cpdef normalize(self)
    cpdef add(self, Position pos)
    cpdef Position max(self)
    cpdef Position sum(self)
    cpdef Position min(self)
    cpdef clear(self)
    cpdef Position overlapWithRange(
        self, Position start, Position end)
    cpdef Position intersectionWithSegments(
        self, SegmentList other, mode = *)
    cpdef void intersect(self, SegmentList other)
    cpdef PositionList clone(self)

    # C only methods
    cdef _resize(self, int npositions)
    cdef _add(self, Position pos)
    cdef int _getInsertionPoint(self, Position other)
    cdef Position overlap(self, Segment other)

    # Functions overloaded from CoordinateList
    cdef off_t toMMAP(self, void *, int, off_t)
    cdef void fromMMAP(self)
