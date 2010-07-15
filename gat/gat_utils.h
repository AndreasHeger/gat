#ifndef __GATUTILS_H
#define __GATUTILS_H

#include <stdlib.h>

// use bisection to return the index of the first element i in base
// such that base[i] <= target. When there is no such index, return
// nmemb.
long searchsorted(void * base,
		  size_t nmemb,
		  size_t size,
		  const void * target,
		  int(*compar)(const void *, const void *));

#endif
		   