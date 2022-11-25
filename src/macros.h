
#ifndef __MACROS_H__
#define __MACROS_H__

#include <string.h>

#ifdef _GNU_SOURCE
#define STRSTRCOMP(A,B) strcasestr((A),(B))
#else
#define STRSTRCOMP(A,B) strstr((A),(B))
#endif

#endif
