# cython: embedsignature=True
# cython: profile=False

cimport cython

from cpython.bytes cimport PyBytes_AsString, PyBytes_FromStringAndSize
from cpython cimport PyBytes_Check, PyUnicode_Check
from libc.stdlib cimport qsort, calloc, malloc, realloc, free
from libc.string cimport memcpy
from libc.errno cimport errno
from posix.types cimport off_t
from posix.mman cimport mmap, munmap, shm_open, shm_unlink
from posix.mman cimport MAP_SHARED, PROT_READ, PROT_WRITE
from posix.stat cimport S_IRUSR, S_IWUSR
from posix.fcntl cimport O_CREAT, O_RDWR, O_RDONLY

from SegmentList cimport Position, Segment, SegmentList

cdef bytes force_bytes(object s, encoding="ascii"):
    """convert string or unicode object to bytes, assuming
    ascii encoding.
    """
    if s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding)
    else:
        raise TypeError("Argument must be string, bytes or unicode.")


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
                 PositionList clone=None,
                 iter=None,
                 sort=False,
                 normalize=False,
                 PositionList share=None,
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
        self.chunk_size = 10000

        # an empty list is normalized and thus sorted
        self.is_normalized = 1
        self.is_sorted = 1
        self.is_shared = False
        self.is_slave = False
        self.shared_fd = -1
        self.key = None
        cdef int idx
        cdef Position pos
        cdef void * retval

        # initialize list from shared memory
        # allocated is 0
        if share != None:
            self.npositions = share.npositions
            if self.npositions > 0:
                assert share.shared_fd >= 0, "sharing from a non-shared SegmentList"

                self.allocated = 0
                self.is_sorted = share.is_sorted
                self.is_normalized = share.is_normalized
                self.is_shared = True
                self.is_slave = True
                self.chunk_size = share.chunk_size
                retval = <Position*>mmap(NULL, 
                                         self.npositions * sizeof(Segment),
                                         PROT_READ, 
                                         MAP_SHARED, 
                                         share.shared_fd, 
                                         0)

                error_number = errno

                if retval == <void*>-1:
                    raise ValueError(
                        "could not read list from shared segment - error=%i" %
                        error_number )
                self.positions = <Position*>retval

        # create from pickled representation
        # if a shared object is pickled, it will be re-constituted 
        # as a slave list.
        elif unreduce:
            self.npositions, self.allocated, self.is_sorted, \
                self.is_normalized, \
                self.chunk_size, self.key = unreduce[:6]
            if type(unreduce[6]) == int:
                # shared memory
                shared_fd = unreduce[6]
                
                retval = mmap(NULL, 
                              self.npositions * sizeof(Position),
                              PROT_READ, 
                              MAP_SHARED, 
                              shared_fd, 
                              0)

                error_number = errno

                if retval == <void*>-1:
                    raise ValueError(
                        "could not unpickle as slave - error=%i" % error_number)

                self.positions = <Position *>retval

                # mark memory as slave - nothing allocated and no shared_fd
                self.allocated = 0
                self.shared_fd = -1
                self.is_shared = True
                self.is_slave = True
            else:
                p = PyBytes_AsString(unreduce[6])
                self.positions = <Position*>malloc(self.npositions * sizeof(Position))
                memcpy(self.positions, p, cython.sizeof(Position) * self.npositions)

        # clone from another list
        elif clone != None:
            self.npositions = self.allocated = clone.npositions
            self.positions = <Position*>calloc(clone.npositions, sizeof(Position))
            memcpy(self.positions,
                   clone.positions,
                   clone.npositions * sizeof(Position))
            self.is_sorted = clone.is_sorted
            self.is_normalized = clone.is_normalized
            
        # create from an iterable
        elif iter:
            a = tuple(iter)
            self.npositions = self.allocated = len(a)
            self.positions = <Position*>calloc(self.npositions, sizeof(Position))
            idx = 0
            for pos in a:
                self.positions[idx] = pos
                idx += 1
            self.is_normalized = False
            self.is_sorted = False

        # allocated memory only, list remains empty
        elif allocate:
            self.allocated = allocate
            self.positions = <Position*>calloc(allocate, sizeof(Position))

        if sort:
            self.sort()

        if normalize:
            self.normalize()

    def __reduce__(self):
        '''pickling function - returns class contents as a tuple.'''

        cdef bytes data

        if self.shared_fd >= 0:
            return (buildPositionList, (self.npositions, 
                                        self.allocated, 
                                        self.is_sorted,
                                        self.is_normalized,
                                        self.chunk_size, 
                                        self.key,
                                        self.shared_fd))

        else:
            data = PyBytes_FromStringAndSize(
                <char*>self.positions, \
                self.npositions * cython.sizeof(Position) * 2)
        
            return (buildPositionList, (self.npositions, 
                                        self.allocated,
                                        self.is_sorted,
                                        self.is_normalized,
                                        self.chunk_size, 
                                        self.key,
                                        data))
    cpdef Position max(self):
        '''return maximum coordinate.'''
        assert self.is_sorted, "PositionList : max from unsorted list"
        if self.npositions == 0:
            return 0
        return self.positions[self.npositions - 1]

    cpdef Position min(self):
        '''return minimum coordinate.'''
        assert self.is_sorted, "PositionList : min from unsorted list"
        if self.npositions == 0:
            return 0
        return self.positions[0]

    cpdef Position sum(self):
        '''return the sum (number) of positions.'''
        return self.npositions

    cpdef clear(self):
        '''removes all positions from list.

        Note that this method does not shrink the allocated memory.
        '''
        self.npositions = 0
        self.is_sorted = True
        self.is_normalized = True

    def __len__(self):
        return self.npositions

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
        self.is_sorted = False
        self.is_normalized = False

    cpdef add(self, Position pos):
        '''add a new position.

        The list will not be normalized automatically - call
        :meth:`PositionList.normalize`.
        '''
        self._add(pos)

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
                self.positions[idx] = start + (end - start) // 2
                idx += 1
        elif method == "start":
            for start, end in segments:
                if start == end:
                    continue
                self.positions[idx] = start
                idx += 1
        elif method == "end":
            for start, end in segments:
                if start == end:
                    continue
                self.positions[idx] = end
                idx += 1
        else:
            raise ValueError("unknow method '%s'" % method)

        if segments.isNormalized:
            self.is_sorted = True
            self.is_normalized = True
        else:
            self.is_sorted = False
            self.is_normalized = False

        self.npositions = idx

    cpdef sort(self):
        """sort points."""
        if self.npositions == 0:
            self.is_sorted = True
            return

        qsort(<void*>self.positions,
              self.npositions,
              sizeof(Position),
              &cmpPosition)

        self.is_sorted = True

    cpdef normalize(self):
        """normalized position list.

        Overlapping positions will be merged.
        """
        if self.npositions == 0:
            self.is_normalized = True
            return

        self.sort()

        cdef int idx = 0
        cdef int insertion_idx = 0
        cdef Position last_position = -1
        cdef Position current_position
        for idx from 0 <= idx < self.npositions:
            current_position = self.positions[idx]
            if current_position != last_position:
                self.positions[insertion_idx] = current_position
                last_position = current_position
                insertion_idx += 1

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
        if other < self.positions[0]:
            return -1
        if other >= self.positions[self.npositions-1]:
            return self.npositions

        return searchsorted(self.positions,
                            self.npositions,
                            sizeof(Position),
                            &other,
                            &cmpPosition)

    cdef Position overlap(self, Segment other):
        """return the number of nucleotides overlapping with other.

        Arguments
        ---------
        other : Segment
            Segment to compute overlap with

        Returns
        ------
        overlap : int
            Number of residues overlapping
        """
        assert self.is_sorted, "PointList is not sorted"

        cdef Position count = 0
        cdef int idx = self._getInsertionPoint(other.start)
        # if beyond last element, check if overlap with
        # previous elements
        if idx == self.npositions:
            idx -= 1
            while idx > 0 and \
                  self.positions[idx] >= other.start and \
                  self.positions[idx] < other.end:
                idx -= 1
                count += 1
            return count

        if idx < 0:
            idx = 0

        while idx < self.npositions and \
              self.positions[idx] >= other.start and \
              self.positions[idx] < other.end:
            count += 1
            idx += 1
        return count
    
    cpdef Position overlapWithRange(self, Position start, Position end):
        """return the number of nucleotides overlapping with a range.

        Arguments
        ---------
        start : int
            Begin of range.
        end : int
            End of range.

        Returns
        ------
        overlap : int
            Number of residues overlapping
        """

        cdef Segment s
        s = Segment(start, end)
        return self.overlap(s)

    cpdef Position intersectionWithSegments(
        self,
        SegmentList other,
        mode="base"):
        """return number of positions overlapping with *other*.

        `mode` is there for compatibility with SegmentList.


        Arguments
        ---------
        other : SegmentList
             SegmentList to overlap with
        mode : string
             Unused.

        Returns
        -------
        overlap : int
            Number of positions
        """

        assert self.is_sorted, "Intersection from unsorted positions"
        assert other.isNormalized, "Intersection with non-normalized segments"

        cdef int this_idx, other_idx, last_this_idx, last_other_idx
        this_idx = other_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Position this_position = -1
        cdef Segment other_segment = Segment(0,0)
        cdef Position noverlap = 0

        while this_idx < self.npositions and other_idx < other.nsegments:
            if last_this_idx != this_idx:
                this_position = self.positions[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # skip segments in this not overlapping other
            if this_position < other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_position:
                other_idx += 1
            else:
                noverlap += 1
                this_idx += 1

        return noverlap

    cpdef void intersect(self, SegmentList other):
        '''intersect this position list with a SegmentList.

        This method removes all Positions that are not contained
        in Segments in `other`.
        '''

        assert self.is_sorted, "Intersection from unsorted positions"
        assert other.isNormalized, "Intersection with non-normalized segments"

        if self.npositions == 0:
            return

        if len(other) == 0:
            self.clear()
            return

        cdef int this_idx, other_idx, last_this_idx, last_other_idx, working_idx
        this_idx = other_idx = working_idx = 0
        last_this_idx = last_other_idx = -1
        cdef Position this_position = -1
        cdef Segment other_segment = Segment(0,0)
        cdef Position noverlap = 0

        while this_idx < self.npositions and other_idx < other.nsegments:
            if last_this_idx != this_idx:
                this_position = self.positions[this_idx]
                last_this_idx = this_idx
            if last_other_idx != other_idx:
                other_segment = other.segments[other_idx]
                last_other_idx = other_idx

            # skip segments in this not overlapping other
            if this_position < other_segment.start:
                this_idx += 1
            # skip segments in other not overlapping this
            elif other_segment.end <= this_position:
                other_idx += 1
            else:
                self.positions[working_idx] = this_position
                working_idx += 1
                this_idx += 1

        self.npositions = working_idx

    def __str__(self):
        return str(self.asList())

    def asList(self):

        cdef int idx
        result = []
        for idx from 0 <= idx < self.npositions:
            result.append(self.positions[idx])
        return result

    cpdef PositionList clone(self):
        """return a copy of self."""
        clone = PositionList(allocate=self.npositions)
        memcpy(clone.positions,
               self.positions,
               self.npositions * sizeof(Position))
        clone.npositions = self.npositions
        clone.is_sorted = self.is_sorted
        clone.is_normalized = self.is_normalized
        return clone

    def __dealloc__(self):
        cdef int fd

        if self.positions != NULL:
            # unshared memory
            if self.is_shared:
                if self.is_slave:
                    pass
                else:
                    # free shared memory as master
                    munmap(self.positions, 
                           self.npositions * sizeof(Position))
            else:
                # free private copy
                free(self.positions)

        if self.shared_fd != -1:
            fd = shm_unlink(self.key)
            if fd == -1:
                raise OSError("could not unlink shared memory")

    cdef void fromMMAP(self):
        '''retrieve data from mmapped location to private copy.
        '''
        # copy data back to unshared memory
        cdef Position * s
        cdef off_t nbytes = sizeof(Position) * self.npositions
    
        if self.npositions == 0:
            return

        s = <Position *>malloc(nbytes)
        if s == NULL:
            raise ValueError( "could not allocate memory when unsharing" )
        
        memcpy(s, self.positions, nbytes)
        
        self.positions = s
        self.allocated = self.npositions
        self.is_shared = False
        self.is_slave = False

    cdef off_t toMMAP(self, 
                      void * mmap, 
                      int shared_fd, 
                      off_t offset):
        '''move contents of list to shared memory 
        at *mmap* with *offset*.
        '''
        if self.shared_fd != -1: 
            raise ValueError("can not move shared PositionList to mmap")

        if self.npositions == 0:
            # do not move empty position lists
            return offset
        
        cdef off_t nbytes = sizeof(Position) * self.npositions
        cdef Position * p = <Position *>mmap + offset

        # copy data
        memcpy(p, self.positions, nbytes)

        # free allocated private memory
        self.allocated = 0
        free(self.positions)
        self.positions = p

        self.is_shared = True
        self.is_slave = True

        return offset + self.npositions

def buildPositionList(*args):
    '''pickling helper function.
    
    Return a re-constructed SegmentList object.
    '''
    return PositionList(unreduce=args)

