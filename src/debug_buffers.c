
#include "debug_buffers.h"
#include "logging.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

/*
  Function to check for NaN or INF in vectors. Vectors must
  have been allocated with NuRe functions to start at index 1.

  float *vector:    pointer to float vector buffer
  int nt:           time step
  size_t len:       length of vector
  int id:           numerical identifier for debug statement
  int min_nt:       start report when nt reaches min_nt; for 
                    nt<min_nt, no debugging takes place
  const char *cbuf: name of vector as char array

  Return:
  bool b_err:       true if at least one NaN or INF has been 
                    found, otherwise false
 */
bool debug_check_vector(float *vector, int nt, size_t len, int id, int min_nt, const char *cbuf) {

  bool b_err = false;

  if (nt<min_nt) { return b_err; }

  int log_level = log_get_level();
  log_set_level(LOG_DEBUG);

  log_debug("(ID=%d): nt=%d, checking vector %s, len %zu\n", id,nt,cbuf,len);
  for (size_t i=1; i<=len; ++i) {
    if (isnan(vector[i]) || isinf(vector[i])) {
      b_err = true;
      log_debug("(ID=%d): nt=%d, i=%zu, %s[i]=%f\n", id,nt,i,cbuf,vector[i]);
    }
  }

  /* revert to previous log level */
  log_set_level(log_level);

  return b_err;
}

/*
  Function to check for NaN or INF in matrices. Matrices must
  have been allocated with NuRe functions to start at index 1.

  float **matrix:   pointer to pointer to float matrix buffer
  int nt:           time step
  size_t lenx:      length of matrix in x-direction
  size_t leny:      length of matrix in y-direction
  int id:           numerical identifier for debug statement
  int min_nt:       start report when nt reaches min_nt; for 
                    nt<min_nt, no debugging takes place
  const char *cbuf: name of matrix as char array

  Return:
  bool b_err:       true if at least one NaN or INF has been 
                    found, otherwise false
 */
bool debug_check_matrix(float **matrix, int nt, size_t lenx, size_t leny, int id, int min_nt, const char *cbuf) {
    
  bool b_err = false;
  
  if (nt<min_nt) { return b_err; }

  int log_level = log_get_level();
  log_set_level(LOG_DEBUG);

  log_debug("(ID=%d): nt=%d, checking matrix %s, lenx %zu, leny %zu\n", id,nt,cbuf,lenx,leny);
  for (size_t i=1; i<=lenx; ++i) {
    for (size_t j=1; j<leny; ++j) {
      if (isnan(matrix[j][i]) || isinf(matrix[j][i])) {
	b_err = true;
	log_debug("(ID=%d): nt=%d, i=%zu, j=%zu, %s[j][i]=%f\n", id,nt,i,j,cbuf,matrix[j][i]);
      }
    }
  }

  /* revert to previous log level */
  log_set_level(log_level);

  return b_err;
}
