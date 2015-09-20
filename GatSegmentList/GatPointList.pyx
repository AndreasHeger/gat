# cython: embedsignature=True
# cython: profile=False

import types, collections, re, os, math, random

from cpython cimport PyString_AsString, PyString_FromStringAndSize
from cpython.version cimport PY_MAJOR_VERSION

cimport cython
from libc.stdint cimport UINT16_MAX


@cython.profile(False)
cdef int cmpPositions(const_void_ptr s1, const_void_ptr s2):
    return <int>(<Position *>s1) - <int>(<Position *>s2)


cdef class SegmentList:
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
        self.segments = NULL
        self.nsegments = 0
        self.allocated = 0

        # an empty list is normalized and thus sorted
        self.is_normalized = 1
        self.is_sorted = 1

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

        self.nsegments = self.allocated = len(segments)
        self.segments = <Point*>calloc(self.nsegments, sizeof(Point))
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

        self.nsegments = idx

    cpdef sort(self):
        """sort points."""
        if self.nsegments == 0:
            return

        qsort(<void*>self.segments,
              self.nsegments,
              sizeof(Segment),
              &cmpPositions)

        self.is_sorted = True

    cpdef normalize(self):
        """normalized Point list.

        Overlapping points will be merged.
        """
        if self.nsegments == 0:
            self.is_normalized = True
            return

        self.sort()

        cdef int idx = 0
        cdef int insertion_idx = 0
        cdef Position last_position
        cdef Position current_position
        for idx from 0 <= idx < self.nsegments
            current_position = self.segments[idx]
            if current_position != last_position:
                self.segments[insertion_idx] = current_position
                last_position = current_position
                insertion_index += 1
            idx += 1

        self.nsegments = insertion_idx
        self._resize(self.nsegments)
        self.is_normalized = True
        
    cdef _resize(self, int nsegments):
        """resize segment list to nsegments."""
        # do not free all if there are no segments
        if nsegments == 0:
            nsegments = 1

        assert nsegments >= self.nsegments, "resizing will loose segments"
        
        if self.allocated == 0:
            self.segments = <Position*>malloc(nsegments * sizeof(Position))
            if not self.segments:
                raise MemoryError("out of memory when allocation %i bytes" %
                                  sizeof(nsegments * sizeof(Position)))
        else:
            self.segments = <Position*>realloc(self.segments,
                                               nsegments * sizeof(Position))
            if not self.segments:
                raise MemoryError("out of memory when allocation %i bytes" %
                                  sizeof(self.nsegments * sizeof(Position)))

        self.allocated = nsegments

