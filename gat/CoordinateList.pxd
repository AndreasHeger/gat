from posix.types cimport off_t

cdef class CoordinateList:

    # Flag indicating if data is residing in shared memory
    cdef bint is_shared
    # Flag if object is owner of data or not.
    cdef bint is_slave
    # Chunk size for allocating an empty object
    cdef int chunk_size
    # File descriptor of shared memory
    cdef int shared_fd
    # Key to identify shared memory segment
    cdef key

    cdef off_t toMMAP(self, void *, int, off_t)
    cdef void fromMMAP(self)
