# cython: embedsignature=True
# cython: profile=False

from libc.stdlib cimport qsort, calloc, malloc, realloc

cimport cython
from libc.stdint cimport UINT16_MAX

from SegmentList cimport Position, SegmentList

# trick to permit const void * in function definitions
cdef extern from *:
    ctypedef void * const_void_ptr "const void*"

@cython.profile(False)
cdef int cmpPosition(const_void_ptr s1, const_void_ptr s2) nogil: 
    return (<Position*>s1)[0] - (<Position*>s2)[0]


cdef class PositionList:
    '''list of segments.

    A list of (non-overlapping) segments.

    The intersection algorithms assume that segments
    are non-overlapping and sorted by start.

    Call normalize to merge overlapping segments.
    '''
    def __init__(self, 
                 int allocate=0,
                 SegmentList clone=None,
                 iter=None,
                 normalize=False,
                 SegmentList share=None,
                 unreduce=None):
        '''create empty list of segments.
        
        Creating is a lazy operation, no memory is allocated unless
        explitely given. A SegmentList might be initialized using
        the following options:

        If *allocate* is set, create an empty list, but allocate slots

        If *clone* is given, create a copy an existing list.

        If *iter* is given, initialize list from iterable.

        If *share* is set, the list will be a slave. It will share
        memory with a _master_, which has been set up by calling the 
        :meth:`share` method. Memory access for the slave is read-only.

        If *unreduce* is given, the list will be re-constituted
        from a pickled data representation (see __reduce__). Note that
        if a _master_ has been pickled, the reconstituted object will be
        a _slave_.

        If *normalize* is set, the list will be normalized.
        '''
        # initialize empty list
        self.positions = NULL
        self.npositions = 0
        self.allocated = 0

        # an empty list is normalized and thus sorted
        self.is_normalized = 1
        self.is_sorted = 1

    cdef _add(self, Position position):
        '''add a new position.

        The list will not be normalized automatically - call
        :meth:`PositionList.normalize`.
        '''

        if self.allocated == 0:
            self.positions = <Position*>malloc(self.chunk_size * sizeof(Position))
            if not self.positions:
                raise MemoryError("out of memory when allocation %i bytes" %
                                  sizeof( self.chunk_size * sizeof(Position)))
            self.allocated = self.chunk_size
        elif self.npositions == self.allocated:
            self.allocated *= 2
            self.positions = <Position*>realloc(self.positions,
                                              self.allocated * sizeof(Position))
            if not self.positions:
                raise MemoryError(
                    "out of memory when allocation %i bytes" %
                    sizeof(self.allocated * sizeof(Position)))

        self.positions[self.npositions] = position
        self.npositions += 1
        self.flag = 0

    def fromSegmentList(self, SegmentList segments, method="midpoint"):
        """fill from SegmentList.

        Valid methods are ``midpoint``, ``start`` and ``end``.

        Arguments
        ---------
        segments : SegmentList
             Segments to take from.
        method : string
             Method to convert segments to points.
        """
        self.clear()

        cdef int idx

        self.npositions = self.allocated = len(segments)
        self.positions = <Position*>calloc(self.npositions, sizeof(Position))
        idx = 0
        if method == "midpoint":
            for start, end in segments:
                if start == end:
                    continue
                self.points[idx] = start + (end - start) // 2
                idx += 1
        elif method == "start":
            for start, end in segments:
                if start == end:
                    continue
                self.points[idx] = start
                idx += 1
        elif method == "end":
            for start, end in segments:
                if start == end:
                    continue
                self.points[idx] = end
                idx += 1

        if segments.isNormalized():
            self.is_sorted = True
            self.is_normalized = True
    
        else:
            self.is_sorted = False
            self.is_normalized = False

        self.npositions = idx

    cpdef sort(self):
        """sort points."""
        if self.npositions == 0:
            return

        qsort(<void*>self.positions,
              self.npositions,
              sizeof(Position),
              &cmpPosition)

        self.is_sorted = True

    cpdef normalize(self):
        """normalized Point list.

        Overlapping points will be merged.
        """
        if self.npositions == 0:
            self.is_normalized = True
            return

        self.sort()

        cdef int idx = 0
        cdef int insertion_idx = 0
        cdef Position last_position
        cdef Position current_position
        for idx from 0 <= idx < self.npositions:
            current_position = self.positions[idx]
            if current_position != last_position:
                self.positions[insertion_idx] = current_position
                last_position = current_position
                insertion_index += 1
            idx += 1

        self.npositions = insertion_idx
        self._resize(self.npositions)
        self.is_normalized = True
        
    cdef _resize(self, int npositions):
        """resize position list to npositions."""
        # do not free all if there are no positions
        if npositions == 0:
            npositions = 1

        assert npositions >= self.npositions, "resizing will loose positions"
        
        if self.allocated == 0:
            self.positions = <Position*>malloc(npositions * sizeof(Position))
            if not self.positions:
                raise MemoryError("out of memory when allocation %i bytes" %
                                  sizeof(npositions * sizeof(Position)))
        else:
            self.positions = <Position*>realloc(self.positions,
                                               npositions * sizeof(Position))
            if not self.positions:
                raise MemoryError("out of memory when allocation %i bytes" %
                                  sizeof(self.npositions * sizeof(Position)))

        self.allocated = npositions

    cdef int _getInsertionPoint(self, Position other):
        """return insertion point.

        The insertion point denotes the element at which
        or after which `other` should be inserted to keep
        the sort order.

        This method does not check if the list is sorted. Peform
        this check before calling _getInsertionPoint.

        Returns
        -------
        position : int
           The function returns -1 if other is before the first element
           or the list is empty. The function returns npositions if other
           is after the last element.
        """

        cdef int idx
        if self.npositions == 0:
            return -1

        # avoid out of range searches
        if other.start >= self.positions[self.npositions-1].end:
            return self.npositions
        if other.end <= self.positions[0].start: 
            return -1

        idx = searchsorted(self.positions,
                           self.npositions,
                           sizeof(Position),
                           &other,
                           &cmpPosition)

        if idx == self.npositions:
            return idx-1
        elif self.positions[idx].start != other.start:
            return idx-1
        else:
            return idx

    cdef Position overlap(self, Position other):
        """return the number of nucleotides overlapping with other.

        Arguments
        ---------
        other : Position
            Position to compute overlap with

        Returns
        ------
        overlap : int
            Number of residues overlapping
        """
        assert self.is_sorted, "PointList is not sorted"

        cdef int idx
        idx = self._getInsertionPoint(other)

        # deal with border cases
        # permit partial overlap
        if idx == self.npositions:
            idx -= 1
        elif idx == -1:
            idx=0

        cdef Position count = 0

        while idx < self.npositions and self.positions[idx] == other:
            count += 1
            idx += 1
        return count
    
