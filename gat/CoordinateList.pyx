# cython: embedsignature=True
# cython: profile=False

cimport cython

from posix.types cimport off_t

cdef class CoordinateList:
    
    def __init__(self):
        self.is_shared = False
        self.is_slave = False

    cdef void fromMMAP(self):
        '''retrieve data from mmapped location to private copy.
        '''
        raise NotImplementedError()

    cdef off_t toMMAP(self, 
                      void * mmap, 
                      int shared_fd, 
                      off_t offset):
        '''move contents of list to shared memory 
        at *mmap* with *offset*.
        '''
        raise NotImplementedError()
