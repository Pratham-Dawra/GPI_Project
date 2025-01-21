
#ifndef DEBUG_BUFFERS_H
#define DEBUG_BUFFERS_H

#include <stdlib.h>
#include <stdbool.h>

bool debug_check_vector(float *, int, size_t, int, int, const char *);
bool debug_check_matrix(float **, int, size_t, size_t, int, int, const char *);

#endif
