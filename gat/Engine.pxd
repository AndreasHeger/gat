from libc.stdio cimport FILE

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

cdef class IntervalContainer:
    cdef int shared_fd     	   
    cdef str shared_fn
    cdef intervals
    cdef str name
 
cdef class IntervalDictionary(IntervalContainer):
    pass

cdef class IntervalCollection(IntervalContainer):
    pass


    
