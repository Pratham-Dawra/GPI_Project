
#ifndef __MACROS_H__
#define __MACROS_H__

#include <math.h>
#include <string.h>

#ifdef _GNU_SOURCE
#define STRSTRCOMP(A,B) strcasestr((A),(B))
#else
#define STRSTRCOMP(A,B) strstr((A),(B))
#endif

#define STRING_SIZE 256

#define NPAR 41       // total max. number of json parameters
#define NSPAR 12      // number of source parameters (matrix size)

#define PI M_PI

#define iround(x) ((int)(floor)((x)+0.5))

#define min(x,y) (((x)<(y))?(x):(y))

#define max(x,y) (((x)<(y))?(y):(x))

#define fsign(x) (((x)<0.0)?(-1):1)

#endif
