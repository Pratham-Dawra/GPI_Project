
#ifndef __UTIL_H_INCLUDED__
#define __UTIL_H_INCLUDED__

/* utility functions */
void dt_mult(int nx, int ny, float dt, float **a);
double maximum(float **a, int nx, int ny);
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

float *vector_c(int nx);
float **matrix_c(int nx, int ny);

#endif
