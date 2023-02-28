
#ifndef __UTIL_H_INCLUDED__
#define __UTIL_H_INCLUDED__

#include <stdlib.h>

/* utility functions */
void dt_mult(int nx, int ny, float dt, float **a);
double maximum(float **a, int nx, int ny);
void shift_var2(float ***var1, float ***var2, float ***var3, float ***var4);
void shift_var3(float ****var1, float ****var2, float ****var3, float ****var4);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

/* =========== Contiguous memory routines in C =========== */

/* To pass a malloc2d'ed buffer to a function like memset or any
 * other function that expects a flat one-dimensional buffer, use
 * FLATBUF2(buffer) as argument for the function. Use FLATBUF3
 * accordingly for malloc3d'ed buffers.
 */
#define FLATBUF2(x) (*(x))      // or use &(x[0][0])
#define FLATBUF3(x) (**(x))     // or use &x[0][0][0])

/* All memory allocated with the following functions can be free'ed with
 * a single call to free(), e.g., buffer = malloc2d(...) and free(buffer). */

/*! Allocation of 1D contiguous memory buffer
 *  @param[in] n1 Length of first dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 */
void *malloc1d(size_t n1, size_t typesize);

/*! Allocation of 2D contiguous memory buffer
 *  @param[in] n1 Length of first dimension
 *  @param[in] n2 Length of second dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 *  @note Address buffer by buffer[i1][i2]; fast axis is n2
 */
void **malloc2d(size_t n1, size_t n2, size_t typesize);

/*! Allocation of 3D contiguous memory buffer
 *  @param[in] n1 Length of first dimension
 *  @param[in] n2 Length of second dimension
 *  @param[in] n3 Length of third dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 *  @note Address buffer by buffer[i1][i2][i3]; fast axis is n3
 */
void ***malloc3d(size_t n1, size_t n2, size_t n3, size_t typesize);

/*! Allocation of 1D contiguous memory buffer, buffer is zeroed
 *  @param[in] n1 Length of first dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 */
void *calloc1d(size_t n1, size_t typesize);

/*! Allocation of 2D contiguous memory buffer, buffer is zeroed
 *  @param[in] n1 Length of first dimension
 *  @param[in] n2 Length of second dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 *  @note Address buffer by buffer[i1][i2]; fast axis is n2
 */
void **calloc2d(size_t n1, size_t n2, size_t typesize);

/*! Allocation of 3D contiguous memory buffer, buffer is zeroed
 *  @param[in] n1 Length of first dimension
 *  @param[in] n2 Length of second dimension
 *  @param[in] n3 Length of third dimension
 *  @param[in] typesize Size of type to allocate, i.e., sizeof(float)
 *  @return Pointer to memory buffer
 *  @note Address buffer by buffer[i1][i2][i3]; fast axis is n3
 */
void ***calloc3d(size_t n1, size_t n2, size_t n3, size_t typesize);

#endif
