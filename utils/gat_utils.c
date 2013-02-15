#include <stdlib.h>
#include <zlib.h>
#include "gat_utils.h"

// return insertion_point, after numpy-1.0.4/numpy/core/src/multiarraymodule.c
// local_search_left
// can be replaced by bsearch?
long searchsorted(const void * base,
		  size_t nmemb,
		  size_t size,
		  const void * target,
		  int(*compar)(const void *, const void *)
		  )
{

  size_t imin = 0;
  size_t imax = nmemb;
  // increment in bytes
  char * b = (char*)base;
  
  while (imin < imax) 
    {
      size_t imid = imin + ((imax - imin) >> 1);
      
      if (compar( &b[imid*size], target) < 0)
	imin = imid + 1;
      else
	imax = imid;
    }

  return imin;
}

// return insertion_point, after numpy-1.0.4/numpy/core/src/multiarraymodule.c
// local_search_left
// can be replaced by bsearch?
long searchargsorted(const void * base,
		     int * sorted,
		     size_t nmemb,
		     size_t size,
		     const void * target,
		     int(*compar)(const void *, const void *)
		     )
{

  size_t imin = 0;
  size_t imax = nmemb;
  // increment in bytes
  char * b = (char*)base;
    
  while (imin < imax) 
    {
      size_t imid = imin + ((imax - imin) >> 1);
      if (compar( &b[sorted[imid]*size], target) < 0)
	imin = imid + 1;
      else
	imax = imid;
    }

  return imin;
}

int toCompressedFile( const unsigned char * buffer, size_t uncompressed_size, FILE * output_f )
{
  uLongf compressed_size = uncompressed_size * 2;

  Bytef * compressed = (Bytef *)calloc( compressed_size, sizeof(Bytef) );

  int level = 9;

  int zok = compress2(compressed, &compressed_size, buffer, uncompressed_size, level);

  if ( zok != Z_OK || fwrite( &compressed_size, sizeof(uLongf), 1, output_f ) != 1 || ferror( output_f ))
    {
      free( compressed );
      return zok;
    }

  if ( fwrite(compressed, 1, compressed_size, output_f) != compressed_size || ferror(output_f))
    {
      free( compressed );
      return -10;
    }
  
  free( compressed );
  return zok;

}


// save compressed data into buffer. Buffer has to be large enough.
int fromCompressedFile( unsigned char * buffer, size_t uncompressed_size, FILE * input_f )
{

  uLongf compressed_size;

  if ( fread( &compressed_size, sizeof(uLongf), 1, input_f) != 1 || ferror(input_f))
    {
      return Z_ERRNO;
    }

  Bytef * compressed = (Bytef *)calloc( compressed_size, sizeof(Bytef) );

  if ( ( fread(compressed, 1, compressed_size, input_f) != compressed_size) || ferror(input_f) )
      {
	free( compressed );
	return Z_ERRNO;
      }
  
  int zok = uncompress (buffer, &uncompressed_size, compressed, compressed_size);
  free(compressed);
  return zok;
}
