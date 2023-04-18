
/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
 * some utility-routines from numerical recipes , see NUMERICAL RECIPES IN C,
 * Press et al., 1990, pp 942
 *
 * ----------------------------------------------------------------------*/

#define NR_END 1
#define FREE_ARG char *

#include "util.h"
#include "logging.h"
#include <math.h>
#include <stdlib.h>

#define UNUSED(x) (void)(x)

void dt_mult(int nx, int ny, float dt, float **a)
{
    int i, j;

    for (i = 1; i <= nx; i++) {
        for (j = 1; j <= ny; j++) {
            a[j][i] = a[j][i] * dt;
        }
    }
}

double maximum(float **a, int nx, int ny)
{
    /* find absolute maximum of array a[1...nx][1...ny] */
    double maxi = 0.0;
    int i, j;

    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            if (fabs((double)a[i][j]) > maxi)
                maxi = fabs((double)a[i][j]);
    return maxi;
}


float average_matrix(float ** matrix, GlobVar *gv)
{
    /* Return average  of Matrix */
    int i,j;
    float local_sum=0;
    float global_sum;
    float buf1=0, buf2=0;
    float average;
    
    for (j=1;j<=gv->NY;j++){
        for (i=1;i<=gv->NX;i++){
            local_sum+=matrix[j][i];
        }
    }
    
    buf1=local_sum;
    buf2=0;
    MPI_Allreduce(&buf1,&buf2, 1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    global_sum=buf2;
    
    average=global_sum/(gv->NXG*gv->NYG);
    
    MPI_Bcast(&average,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    
    return average;
}

void shift_var2(float ***var1, float ***var2, float ***var3, float ***var4)
{
    /* Adam Bashforth time shift */
    float **shift_1 = *var4;

    *var4 = *var3;
    *var3 = *var2;
    *var2 = *var1;
    *var1 = shift_1;
}

void shift_var3(float ****var1, float ****var2, float ****var3, float ****var4)
{
    /* Adam Bashforth time shift */
    float ***shift_1 = *var4;

    *var4 = *var3;
    *var3 = *var2;
    *var2 = *var1;
    *var1 = shift_1;
}

float *vector(int nl, int nh)
{
    /* allocate a float vector with subscript range v[nl..nh] and initializing
     * this vector, eg. vector[nl..nh]=0.0 */
    float *v;
    int i;

    v = (float *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(float)));
    if (!v)
        log_fatal("Allocation failure in function vector().\n");
    for (i = 0; i < (nh - nl + 1 + NR_END); i++)
        v[i] = 0.0;
    return v - nl + NR_END;
}

int *ivector(int nl, int nh)
{
    /* allocate an int vector with subscript range v[nl..nh] and initializing
     * this vector, eg. ivector[nl..nh]=0 */
    int *v;
    int i;

    v = (int *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v)
        log_fatal("Allocation failure in function ivector().\n");
    for (i = 0; i < (nh - nl + 1 + NR_END); i++)
        v[i] = 0;
    return v - nl + NR_END;
}

unsigned short int *usvector(int nl, int nh)
{
    /* allocate an short int vector with subscript range v[nl..nh] and initializing
     * this vector, eg. ivector[nl..nh]=0 */
    unsigned short int *v;
    int i;

    v = (unsigned short int *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(unsigned short int)));
    if (!v)
        log_fatal("Allocation failure in function usvector().\n");
    for (i = 0; i < (nh - nl + 1 + NR_END); i++)
        v[i] = 0;
    return v - nl + NR_END;
}

unsigned char *cvector(int nl, int nh)
{
    /* allocate an unsigned char vector with subscript range v[nl..nh] */
    unsigned char *v;

    v = (unsigned char *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
    if (!v)
        log_fatal("Allocation failure in function cvector().\n");
    return v - nl + NR_END;
}

unsigned long *lvector(int nl, int nh)
{
    /* allocate an unsigned long vector with subscript range v[nl..nh] and
     * initializing this vector, eg. vector[nl..nh]=0.0 */
    unsigned long *v;
    int i;

    v = (unsigned long *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(unsigned long)));
    if (!v)
        log_fatal("Allocation failure in function lvector().\n");
    for (i = 0; i < (nh - nl + 1 + NR_END); i++)
        v[i] = 0;
    return v - nl + NR_END;
}

double *dvector(int nl, int nh)
{
    /* allocate a double vector with subscript range v[nl..nh] and initializing
     * this vector, eg. vector[nl..nh]=0.0 */

    double *v;
    int i;

    v = (double *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        log_fatal("Allocation failure in function dvector().\n");
    for (i = 0; i < (nh - nl + 1 + NR_END); i++)
        v[i] = 0.0;
    return v - nl + NR_END;
}

float **fmatrix(int nrl, int nrh, int ncl, int nch)
{
    /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    float **m;

    /* allocate pointers to rows */
    m = (float **)malloc((size_t) ((nrow + NR_END) * sizeof(float *)));
    if (!m)
        log_fatal("Allocation failure 1 in function fmatrix().\n");
    m += NR_END;
    m -= nrl;

    /* allocation rows and set pointers to them */
    m[nrl] = (float *)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(float)));
    if (!m[nrl])
        log_fatal("Allocation failure 2 in function fmatrix().\n");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* initializing matrix */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            m[i][j] = 0.0f;

    /* return pointer to array of pointer to rows */
    return m;
}

float **matrix(int nrl, int nrh, int ncl, int nch)
{
    /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    float **m;

    /* allocate pointers to rows */
    m = (float **)malloc((size_t) ((nrow + NR_END) * sizeof(float *)));
    if (!m)
        log_fatal("Allocation failure 1 in function matrix().\n");
    m += NR_END;
    m -= nrl;

    /* allocation rows and set pointers to them */
    m[nrl] = (float *)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(float)));
    if (!m[nrl])
        log_fatal("Allocation failure 2 in function matrix().\n");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* initializing matrix */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            m[i][j] = 0.0f;

    /* return pointer to array of pointer to rows */
    return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m = (double **)malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        log_fatal("Allocation failure 1 in function matrix().\n");
    m += NR_END;
    m -= nrl;

    /* allocation rows and set pointers to them */
    m[nrl] = (double *)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        log_fatal("Allocation failure 2 in function dmatrix().\n");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* initializing matrix */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            m[i][j] = 0.0;

    /* return pointer to array of pointer to rows */
    return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    /* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int **m;

    /* allocate pointers to rows */
    m = (int **)malloc((size_t) ((nrow + NR_END) * sizeof(int *)));
    if (!m)
        log_fatal("Allocation failure 1 in function imatrix().\n");
    m += NR_END;
    m -= nrl;

    /* allocation rows and set pointers to them */
    m[nrl] = (int *)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int)));
    if (!m[nrl])
        log_fatal("Allocation failure 2 in function imatrix().\n");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* initializing matrix */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            m[i][j] = 0;

    /* return pointer to array of pointer to rows */
    return m;
}

unsigned short int **usmatrix(int nrl, int nrh, int ncl, int nch)
{
    /* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
    int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    unsigned short int **m;

    /* allocate pointers to rows */
    m = (unsigned short int **)malloc((size_t) ((nrow + NR_END) * sizeof(unsigned short int *)));
    if (!m)
        log_fatal("Allocation failure 1 in function usmatrix().\n");
    m += NR_END;
    m -= nrl;

    /* allocation rows and set pointers to them */
    m[nrl] = (unsigned short int *)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(unsigned short int)));
    if (!m[nrl])
        log_fatal("Allocation failure 2 in function usmatrix().\n");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* initializing matrix */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            m[i][j] = 0;

    /* return pointer to array of pointer to rows */
    return m;
}

float ***f3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    /* allocate a float 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh]=0.0 */
    int i, j, d, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    float ***t = NULL;

    /* allocate pointers to pointers to rows */
    t = (float ***)malloc((size_t) ((nrow + NR_END) * sizeof(float **)));
    if (!t)
        log_fatal("Allocation failure 1 in function f3tensor().\n");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (float **)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(float *)));
    if (!t[nrl])
        log_fatal("Allocation failure 2 in function f3tensor().\n");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (float *)malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(float)));
    if (!t[nrl][ncl])
        log_fatal("Allocation failure 3 in function f3tensor().\n");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* initializing 3tensor */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            for (d = ndl; d <= ndh; d++)
                t[i][j][d] = 0.0f;

    /* return pointer to array of pointer to rows */
    return t;
}

int ***i3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    /* allocate a integer 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh]=0.0 */
    int i, j, d, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    int ***t = NULL;

    /* allocate pointers to pointers to rows */
    t = (int ***)malloc((size_t) ((nrow + NR_END) * sizeof(int **)));
    if (!t)
        log_fatal("Allocation failure 1 in function i3tensor().\n");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (int **)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int *)));
    if (!t[nrl])
        log_fatal("Allocation failure 2 in function i3tensor().\n");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (int *)malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(int)));
    if (!t[nrl][ncl])
        log_fatal("Allocation failure 3 in function i3tensor().\n");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* initializing 3tensor */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            for (d = ndl; d <= ndh; d++)
                t[i][j][d] = 0;

    /* return pointer to array of pointer to rows */
    return t;
}

float ****f4tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh, int nvl, int nvh)
{
    /* allocate a float 3tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]
     * and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]=0.0 */
    int i, j, d, v, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1, nval = nvh - nvl + 1;
    float ****t = NULL;

    /* allocate pointers to pointers to rows */
    t = (float ****)malloc((size_t) ((nrow + NR_END) * sizeof(float **)));
    if (!t)
        log_fatal("Allocation failure 1 in function f4tensor().\n");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (float ***)malloc((size_t) ((nrow * ncol + NR_END) * sizeof(float *)));
    if (!t[nrl])
        log_fatal("Allocation failure 2 in function f4tensor().\n");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (float **)malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(float)));
    if (!t[nrl][ncl])
        log_fatal("Allocation failure 3 in function f4tensor().\n");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    /* allocate values and set pointers to them */
    t[nrl][ncl][ndl] = (float *)malloc((size_t) ((nrow * ncol * ndep * nval + NR_END) * sizeof(float)));
    if (!t[nrl][ncl][ndl])
        log_fatal("Allocation failure 4 in function f4tensor().\n");
    t[nrl][ncl][ndl] += NR_END;
    t[nrl][ncl][ndl] -= nvl;

    /*for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep; */
    for (d = ndl + 1; d <= ndh; d++)
        t[nrl][ncl][d] = t[nrl][ncl][d - 1] + nval;

    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol * ndep * nval;
        for (j = ncl + 1; j <= nch; j++) {
            t[i][j] = t[i][j - 1] + ndep;
            t[i][j][ndl] = t[i][j - 1][ndl] + ndep * nval;
            for (d = ndl + 1; d <= ndh; d++) {
                t[i][j][d] = t[i][j][d - 1] + nval;
            }
        }
    }

    /* initializing 4tensor */
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
            for (d = ndl; d <= ndh; d++)
                for (v = nvl; v <= nvh; v++)
                    t[i][j][d][v] = 0.0f;

    /* return pointer to array of pointer to rows */
    return t;
}

void free_vector(float *v, int nl, int nh)
{
    /* free a float vector allocated with vector() */
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED(nh);
}

void free_ivector(int *v, int nl, int nh)
{
    /* free a int vector allocated with vector() */
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED(nh);
}

void free_cvector(char *v, int nl, int nh)
{
    /* free a char vector allocated with vector() */
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED(nh);
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
    /* free a float matrix allocated by matrix() */
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED(nrh);
    UNUSED(nch);
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
    /* free a integer matrix allocated by imatrix() */
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED(nrh);
    UNUSED(nch);
}

void free_usmatrix(unsigned short int **m, int nrl, int nrh, int ncl, int nch)
{
    /* free a integer matrix allocated by imatrix() */
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED(nrh);
    UNUSED(nch);
}

void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    /* free a float matrix allocated by f3tensor() */
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
    UNUSED(nrh);
    UNUSED(nch);
    UNUSED(ndh);
}

void free_i3tensor(int ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    /* free a float matrix allocated by i3tensor() */
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
    UNUSED(nrh);
    UNUSED(nch);
    UNUSED(ndh);
}

/* =========== Contiguous memory routines =========== */

void *malloc1d(size_t n1, size_t typesize)
{
    /* allocate memory for vector[0...n1-1] */

    size_t len = n1 * typesize;

    void *buffer = malloc(len);
    if (!buffer)
        log_fatal("malloc of size %zu failed in malloc1d().\n", len);

    return buffer;
}

void *calloc1d(size_t n1, size_t typesize)
{
    void *buffer = calloc(n1, typesize);
    if (!buffer)
        log_fatal("calloc of size %zu failed in calloc1d().\n", n1 * typesize);

    return buffer;
}

void **malloc2d(size_t n1, size_t n2, size_t typesize)
{
    /* allocate memory for matrix[0...n1-1][0...n2-1] with the fast axis being n2;
     * allocation is handled as single malloc to avoid memory fragmentation */

    size_t stride = n2 * typesize;
    void **buffer = malloc(n1 * sizeof(void *) + n1 * stride);
    if (!buffer)
        log_fatal("malloc of size %zu failed in malloc2d().\n", n1 * sizeof(void *) + n1 * stride);

    unsigned char *ptr = (unsigned char *)(buffer + n1);
    for (size_t i = 0; i < n1; i++) {
        buffer[i] = &ptr[i * stride];
    }

    return buffer;
}

void **calloc2d(size_t n1, size_t n2, size_t typesize)
{
    size_t stride = n2 * typesize;

    void **buffer = calloc(n1, sizeof(void *) + stride);
    if (!buffer)
        log_fatal("calloc of size %zu failed in calloc2d().\n", n1 * (sizeof(void *) + stride));

    unsigned char *ptr = (unsigned char *)(buffer + n1);
    for (size_t i = 0; i < n1; i++) {
        buffer[i] = &ptr[i * stride];
    }

    return buffer;
}

void ***malloc3d(size_t n1, size_t n2, size_t n3, size_t typesize)
{
    /* allocate memory for cube[0...n1-1][0...n2-1][0...n3-1] with the fast axis 
     * being n3; allocation is handled as single malloc to avoid memory fragmentation */

    size_t stride1 = n2 * n3 * typesize;
    size_t stride2 = n3 * typesize;
    void ***buffer = malloc(n1 * sizeof(void **) + n1 * n2 * sizeof(void *) + n1 * stride1);
    if (!buffer)
        log_fatal("malloc of size %zu failed in malloc3d().\n", n1 * (sizeof(void **) + n2 * sizeof(void *) + stride1));

    void **p2 = (void **)(buffer + n1);
    unsigned char *p3 = (unsigned char *)(p2 + n1 * n2);
    for (size_t i = 0; i < n1; i++) {
        buffer[i] = &p2[i * n2];
    }
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
            p2[i * n2 + j] = &p3[i * stride1 + j * stride2];
        }
    }

    return buffer;
}

void ***calloc3d(size_t n1, size_t n2, size_t n3, size_t typesize)
{
    size_t stride1 = n2 * n3 * typesize;
    size_t stride2 = n3 * typesize;
    void ***buffer = calloc1d(n1, sizeof(void **) + n2 * sizeof(void *) + stride1);
    if (!buffer)
        log_fatal("calloc of size %zu failed in calloc3d().\n", n1 * (sizeof(void **) + n2 * sizeof(void *) + stride1));

    void **p2 = (void **)(buffer + n1);
    unsigned char *p3 = (unsigned char *)(p2 + n1 * n2);
    for (size_t i = 0; i < n1; i++) {
        buffer[i] = &p2[i * n2];
    }
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
            p2[i * n2 + j] = &p3[i * stride1 + j * stride2];
        }
    }

    return buffer;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef UTIL_C_MAIN

#include <assert.h>

// test main program
int main(void)
{
    size_t n1 = 7;
    size_t n2 = 5;
    size_t n3 = 3;

    float *data1 = (float *)malloc1d(n1, sizeof(float));
    assert(data1);
    size_t count = 0;
    for (size_t i = 0; i < n1; ++i) {
        data1[i] = count++;
    }
    for (size_t i = 0; i < n1; ++i) {
        printf("n1: %zd, address: %p, value: %f\n", i, &(data1[i]), data1[i]);
    }
    free(data1);

    printf("==============================\n");

    data1 = (float *)calloc1d(n1, sizeof(float));
    for (size_t i = 0; i < n1; ++i) {
        printf("n1: %zd, address: %p, value: %f\n", i, &(data1[i]), data1[i]);
    }
    free(data1);

    printf("==============================\n");

    float **data2 = (float **)malloc2d(n1, n2, sizeof(float));
    assert(data2);
    count = 0;
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            data2[i][j] = count++;  // OR *(*(buffer+i)+j) = count++
        }
    }
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            printf("n1: %zd, n2: %zd, address: %p, value: %f\n", i, j, &(data2[i][j]), data2[i][j]);
        }
    }
    free(data2);

    printf("==============================\n");

    data2 = (float **)calloc2d(n1, n2, sizeof(float));
    assert(data2);
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            printf("n1: %zd, n2: %zd, address: %p, value: %f\n", i, j, &(data2[i][j]), data2[i][j]);
        }
    }
    free(data2);

    printf("==============================\n");

    float ***data3 = (float ***)malloc3d(n1, n2, n3, sizeof(float));
    assert(data3);
    count = 0;
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                data3[i][j][k] = count++;   // OR *(*(*(buffer+i)+j)+k) = count++
            }
        }
    }
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                printf("n1: %ld, n2: %ld, n3: %ld, address: %p, value: %f\n", i, j, k, &(data3[i][j][k]),
                       data3[i][j][k]);
            }
        }
    }
    free(data3);

    printf("==============================\n");

    data3 = (float ***)calloc3d(n1, n2, n3, sizeof(float));
    assert(data3);

    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t k = 0; k < n3; ++k) {
                printf("n1: %ld, n2: %ld, n3: %ld, address: %p, value: %f\n", i, j, k, &(data3[i][j][k]),
                       data3[i][j][k]);
            }
        }
    }
    free(data3);

    return 0;
}

#endif

///////////////////////////////////////////////////////////////////////////////
