
#ifndef __READ_SRCSIG_H__
#define __READ_SRCSIG_H__

#include "globvar_struct.h"

/*! Read source signature file.
 *  @param[out] ns number of samples
 *  @param[in] gv global variable structure
 *  @return float buffer with source signature values
 */
float *read_srcsig(int *ns, const GlobVar * gv);

#endif
