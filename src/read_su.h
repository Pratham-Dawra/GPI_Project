
#ifndef __READ_SU_H__
#define __READ_SU_H__

#include "su_gather.h"
#include "su_struct.h"
#include <stdio.h>
#include <stdbool.h>

/* Get number of traces in SU file.
 * @param[in] filep File pointer to SU file (previously opened)
 * @param[out] ns Number of samples per trace
 * @param[out] dt Sampling interval in microseconds
 * @param[out] delrt Delay recording time in milliseconds
 * @return number of traces in file
 * @note File pointer is at beginning of file after calling this function.
 */
size_t su_get_nt(FILE *filep, unsigned short *ns, unsigned short *dt, short *delrt);

/* Read a single seismic trace from disk.
 * @param[in] filep File pointer to SU file (previously opened)
 * @param[in] n Trace number to read in range [0,nt) - zero-based(!)
 * @param[in] ns Number of samples per trace
 * @param[in] b_seek If true, seek to correct file position, otherwise assume file pointer is correct
 * @param[in,out] header SU trace header; if NULL on input, no header is read
 * @param[out] data trace data
 * @return Number of traces read (1)
 */
int su_read_trace(FILE *filep, size_t n, unsigned short ns, bool b_seek, SUhead *header, float *data);

/* Read entire SU file from disk; memory for gather is allocated internally.
 * @param[in] filep File pointer to SU file (previously opened)
 * @param[out] gather struct of type SUgather containing trace headers (nt) and data (nt,ns)
 * @return Number of traces read
 * @note The client is responsible to deallocate memory; call free_SUgather(gather);
 */
size_t su_read_file(FILE *filep, SUgather *gather);

#endif
