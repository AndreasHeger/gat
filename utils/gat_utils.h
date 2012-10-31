#ifndef __GATUTILS_H
#define __GATUTILS_H

#include <stdlib.h>
#include <stdio.h>

// use bisection to return the index of the first element i in base
// such that base[i] <= target. When there is no such index, return
// nmemb.
long searchsorted(const void * base,
		  size_t nmemb,
		  size_t size,
		  const void * target,
		  int(*compar)(const void *, const void *));

long searchargsorted(const void * base,
		     int * sorted,
		     size_t nmemb,
		     size_t size,
		     const void * target,
		     int(*compar)(const void *, const void *));

int toCompressedFile( const unsigned char *, size_t, FILE *);

int fromCompressedFile( unsigned char *, size_t, FILE *);

#endif
		   
