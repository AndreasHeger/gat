#include <stdlib.h>
#include "gat_utils.h"

// return insertion_point, after numpy-1.0.4/numpy/core/src/multiarraymodule.c
// local_search_left
long searchsorted(void * base,
		  size_t nmemb,
		  size_t size,
		  const void * target,
		  int(*compar)(const void *, const void *)
		  )
{

  size_t imin = 0;
  size_t imax = nmemb;
  
  while (imin < imax) 
    {
      size_t imid = imin + ((imax - imin) >> 2);
      
      if (compar( &base[imid*size], target) < 0)
	imin = imid + 1;
      else
	imax = imid;
    }

  return imin;
}

