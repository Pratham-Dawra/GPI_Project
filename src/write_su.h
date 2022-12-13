
#ifndef __WRITE_SU_H__
#define __WRITE_SU_H__

#include "su_gather.h"
#include "su_struct.h"
#include <stdio.h>
#include <stdbool.h>

/*! Write a single seismic trace to disk.
 *  @param[in] filep File pointer to SU file (previously opened)
 *  @param[in] header SU trace header
 *  @param[in] data trace data
 *  @return 0 on success, otherwise a value >0
 *  @note The function logs errors but does not abort. ns must be set correctly in header.
 */
int su_write_trace(FILE * filep, const SUhead * header, const float *data);

/*! Write entire SUgather to disk.
 *  @param[in] filep File pointer to SU file (previously opened).
 *  @param[in] gather struct of type SUgather containing trace headers (nt) and data (nt,ns)
 *  @return 0 on success, otherwise a value >0
 *  @note The function logs errors but does not abort.
 */
int su_write_gather(FILE * filep, const SUgather * gather);

#endif
