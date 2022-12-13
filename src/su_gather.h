
#ifndef __SU_GATHER_H__
#define __SU_GATHER_H__

#include "su_struct.h"
#include <stdlib.h>

typedef struct {
    size_t nt;                  //<! Number of traces in gather
    unsigned short ns;          //<! Number of samples per trace
    SUhead *header;             //<! Trace headers [0..nt)
    float **data;               //<! Traces [0..nt)[0..ns)
} SUgather;

/*! Allocate memory for a seismic gather.
 *  @param[in,out] gather struct of type SUgather
 *  @param[in] nt Number of traces
 *  @param[in] ns Number of samples per trace
 */
void malloc_SUgather(SUgather *gather, size_t nt, unsigned short ns);

/*! Free memory previously allocated for a seismic gather.
 *  @param[in,out] gather struct of type SUgather
 */
void free_SUgather(SUgather *gather);

#endif
