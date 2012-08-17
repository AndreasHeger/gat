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
cdef class SegmentList:

    cdef Segment * segments
    cdef size_t nsegments
    cdef size_t allocated
    cdef int is_normalized
    cdef int chunk_size

    cpdef sort( self )
    cpdef SegmentList extend( self, SegmentList other )
    cpdef add( self, Position start, Position end )
    cpdef trim_ends( self, Position pos, Position size, int forward  )
    cpdef trim( self, Position pos, Position size )
    cpdef normalize( self )
    cpdef merge( self, PositionDifference distance  )
    cpdef check( self )
    cpdef Position getRandomPosition( self )
    cpdef Position overlapWithRange( self, Position start, Position end )
    cpdef SegmentList getOverlappingSegmentsWithRange( self, Position start, Position end )
    cpdef Position overlapWithSegments( self, SegmentList other )
    cpdef Position intersectionWithSegments( self, SegmentList other, mode = * )
    cpdef SegmentList getFilledSegmentsFromStart( self, Position start, PositionDifference remainder )
    cpdef SegmentList getFilledSegmentsFromEnd( self, Position end, PositionDifference remainder )
    cpdef SegmentList filter( self, SegmentList other )
    cpdef SegmentList intersect( self, SegmentList other )
    cpdef extend_segments( self, int extension )
    cpdef expand_segments( self, double expansion )
    cpdef Position sum( self )
    cpdef Position max( self )
    cpdef Position min( self )
    cpdef clear( self )

    cdef _add( self, Segment segment )
    cdef _resize( self, int nsegments )
    cdef insert( self, int idx, Segment seg )
    cdef int _getInsertionPoint( self, Segment other )
    cdef Position overlap( self, Segment other )
    cdef SegmentList getOverlappingSegments( self, Segment other )
    cdef SegmentList truncate( self, Segment other )
