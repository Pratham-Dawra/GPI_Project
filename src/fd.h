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
/* $Id: fd.h 865 2015-09-22 12:57:11Z tmetz $ */
/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD program sofi2D
 *
 *  ---------------------------------------------------------------------*/

#ifndef FD_H_INCLUDED
#define FD_H_INCLUDED

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "globvar_struct.h"

#define iround(x) ((int)(floor)((x)+0.5))
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)<(y))?(y):(x))
#define fsign(x) (((x)<0.0)?(-1):1)

#define PI (3.141592653589793238462643383279502884197169)
#define NPAR 41
#define NSPAR 12     // number of source parameters (matrix size)
//#define STRING_SIZE 74 // previous value, sometimes not enough to handle longer file names
//#define STRING_SIZE 256
#define REQUEST_COUNT 4


/* declaration of functions */
void abs_update_s(int i, int j, float **sxx,float **sxy,float **syy, float **absorb_coeff);

void abs_update_v(int i, int j, float **vx,float **vy, float **absorb_coeff);

void absorb(float **absorb_coeff, GlobVar *gv);

void av_mat(float   **pi, float   **u,
            float   **ppijm, float   **puip, float **pujm);

void av_mue(float **u, float **uipjp, GlobVar *gv);

void av_rho(float **rho, float **rip, float **rjp, GlobVar *gv);

void av_tau(float **taus, float **tausipjp, GlobVar *gv);

void check_fs(GlobVar *gv);

void checkfd(float **prho, float **ppi, float **pu,
             float **ptaus, float **ptaup, float *peta, float *hc, float **srcpos, int nsrc, int **recpos, int ntr, GlobVar *gv);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);

void cpml_update_s_x(int i, int j,float   *vxx, float *vyx,float *K_x, float *a_x,
                     float *b_x, float *K_x_half, float *a_x_half, float *b_x_half ,float **psi_vxx,float **psi_vyx);

void cpml_update_s_y(int i, int j,float *vxy,float *vyy,float *K_y, float *a_y,
                     float *b_y, float *K_y_half, float *a_y_half, float *b_y_half ,float **psi_vyy,float **psi_vxy);

void cpml_update_v_x(int i, int j,float   *sxx_x, float *sxy_x,float *K_x, float *a_x,float *b_x,
                     float *K_x_half, float *a_x_half, float *b_x_half ,float **psi_sxx_x,float **psi_sxy_x);

void cpml_update_v_y(int i, int j,float *sxy_y,float *syy_y,float *K_y, float *a_y,
                     float *b_y, float *K_y_half, float *a_y_half, float *b_y_half ,float **psi_syy_y,float **psi_sxy_y);

void exchange_v(int nd, float **vx, float **vy,
                float **bufferlef_to_rig, float **bufferrig_to_lef,
                float **buffertop_to_bot, float **bufferbot_to_top, GlobVar *gv);

void exchange_s(int nd, float **sxx, float **syy,
                float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
                float **buffertop_to_bot, float **bufferbot_to_top, GlobVar *gv);

void exchange_par(GlobVar *gv);

float *holbergcoeff(GlobVar *gv);

void initfd(GlobVar *gv);

void initproc(GlobVar *gv);

void model_visco(float    **rho, float   **pi, float   **u,
                 float   **taus, float   **taup, float   *eta, GlobVar *gv);

void model_elastic(float    **rho, float   **pi, float   **u, GlobVar *gv);

void model_elastic_VTI(float    **rho, float   **pc11, float   **pc33, float   **pc13, float   **pc55, GlobVar *gv);

void model_elastic_TTI(float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55, float **  pc15,
                       float **  pc35, GlobVar *gv);

void model_visco_vti(float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55,
                     float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55, float *  eta, GlobVar *gv);

void model_visco_tti(float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55, float **  pc15, float **  pc35,
                     float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55, float ** ptau15, float ** ptau35, float *  eta, GlobVar *gv);


void model_ani(float    **rho, float   **c11, float   **c15, float   **c13,
               float   **c35, float   **c33, float   **c55,
               float   **taus, float   **taup, float   *eta);

void matcopy(float **prho, float **ppi, float **pu, float **ptaup,
             float **ptaus, GlobVar *gv);

void matcopy_elastic(float **prho, float **ppi, float **pu, GlobVar *gv);

void matcopy_ani(float **rho, float   **c11, float   **c15, float   **c13,
                 float   **c35, float   **c33, float   **c55, float **taus,
                 float **taup);

void merge(int nsnap, int type, GlobVar *gv);

void mergemod(const char* modfile, int format, GlobVar *gv);

void  outseis_glob(FILE *fpdata, float **section,
                   int **recpos, int ntr, float **srcpos_loc,
                   int ns, int seis_form, int ishot, int comp, GlobVar *gv);

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void operator_s_fd2(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_s_fd4(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_s_fd6(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_s_fd8(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_s_fd10(int i, int j,float   *vxx, float *vyx,float *vxy,
                     float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_s_fd12(int i, int j,float   *vxx, float *vyx,float *vxy,
                     float *vyy, float **vx, float **vy,float *hc, GlobVar *gv);

void operator_v_fd2(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd4(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd6(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd8(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd10(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                     float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd12(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                     float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void par_mult_dt(float **pi, float **u, float **uipjp);

void prepare_update_s(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                      float **puipjp, float **ppi, float **ptaus, float **ptaup,
                      float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                      float *cip, float *cjm, float ***dip, float ***d, float ***e, GlobVar *gv);

void prepare_update_s_4(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                        float **puipjp, float **ppi, float **ptaus, float **ptaup,
                        float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                        float *cip, float *cjm, float ***dip, float ***d, float ***e, GlobVar *gv);

void prepare_update_s_vti(float *peta, float ** pc11, float **pc13, float ** pc33, float **pc55ipjp,
                          float **ptau11, float **ptau13, float ** ptau33, float **ptau55ipjp,
                          float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                          float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                          float *bip, float *cip, GlobVar *gv);

void prepare_update_s_tti(float *peta,
                          float ** pc11, float **pc33, float ** pc13, float **pc55, float **pc15, float **pc35,
                          float ** pc55ipjp, float ** pc15ipjp, float ** pc35ipjp,
                          float **ptau11, float **ptau33, float ** ptau13, float ** ptau55, float ** ptau15, float ** ptau35,
                          float **ptau55ipjp, float **ptau15ipjp,float **ptau35ipjp,
                           float ** pc11u, float **pc33u, float **pc13u, float ** pc55u, float ** pc15u, float ** pc35u,
                          float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                          float *** pc11d, float ***pc33d, float ***pc13d, float *** pc55d,
                          float *** pc15d, float *** pc35d,
                          float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                          float *bip, float *cip, GlobVar *gv) ;

void psource(int nt, float **sxx, float **syy,
             float   **srcpos_loc, float **signals, int nsrc, GlobVar *gv);

void psource_rsg(int nt, float **sxx, float **syy,
                 float   **srcpos_loc, float **signals, int nsrc, GlobVar *gv);

float readdsk(FILE *fp_in, int format);

void readbufs(float **sxx, float **syy,
              float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
              float **buffertop_to_bot, float **bufferbot_to_top);

void readbufv(float **vx, float **vy,
              float **bufferlef_to_rig, float **bufferrig_to_lef,
              float **buffertop_to_bot, float **bufferbot_to_top);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float   **vx, float **vy, float **sxx, float **syy, float **sxy, GlobVar *gv);

void read_par_json(const char *fileinp, GlobVar *gv);


void readmod_visco(float    **rho, float   **pi, float   **u,
                   float   **taus, float   **taup, float   *eta, GlobVar *gv);

void readmod_elastic(float    **rho, float   **pi, float   **u, GlobVar *gv);

void readmod_elastic_vti (float  **  rho, float **  pc11, float **  pc33, float **  pc13,float **  pc55, GlobVar *gv);

void readmod_visco_vti (float  **  rho,
                        float **  pc11, float **  pc33, float **  pc13,float **  pc55,
                        float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55,
                        float *  eta, GlobVar *gv);

void readmod_visco_tti (float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55,
                        float **  pc15, float **  pc35,
                        float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55,
                        float ** ptau15, float ** ptau35, float *  eta, GlobVar *gv);


void readmod_elastic_tti (float  **  rho, float **  pc11, float **  pc33, float **  pc13,float **  pc55, float **  pc15,float **  pc35, GlobVar *gv );

int **receiver(int *ntr, GlobVar *gv);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float   **vx, float **vy, float **sxx, float **syy, float **sxy, GlobVar *gv);

void saveseis_glob(float **sectiondata, int  **recpos, 
                   int ntr, float **srcpos, int ishot,int ns, int sectiondatatype, GlobVar *gv);

void snap(int nt, int nsnap, float **vx, float **vy, float **sxx,
          float **syy, float **u, float **pi, float *hc, GlobVar *gv);

void snapmerge(int nsnap);

float **sources(int *nsrc, GlobVar *gv);

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx,
            float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
            float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu,
                float *hc, GlobVar *gv);

void seismo_rsg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch, GlobVar *gv);

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc, GlobVar *gv);

void subgrid_bounds(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, GlobVar *gv);

void surface(int ndepth, float **pvx, float **pvy,
             float **psxx, float **psyy,
             float **psxy, float *** pp, float *** pq,
             float    **ppi, float    **pu, float **ptaup,
             float **ptaus, float *etajm, float *peta, float *hc, float *K_x, float *a_x, float *b_x, float **psi_vxx, GlobVar *gv);

void surface_elastic(int ndepth, float **vx, float **vy, float **sxx, float **syy,
                     float **sxy, float    **pi, float    **u, float *hc, float *K_x, float *a_x, float *b_x, float**
                     psi_vxx, GlobVar *gv);

void update_s_visc_abs(int *gx, int *gy,
                       float **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float ***r, float *** p, float ***q,
                       float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                       float *cjm, float ***d, float ***e, float ***dip,
                       float **absorb_coeff, GlobVar *gv);

void update_s_visc_abs_4(int *gx, int *gy, int nt,
                         float **vx, float **vy, float **sxx, float **syy,
                         float **sxy, float ***r, float *** p, float ***q,
			 float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                         float *cjm, float ***d, float ***e, float ***dip,
                         float **absorb_coeff, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,
			 float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,
			 float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float ***pr_2,float ***pr_3,float ***pr_4, 
			 float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv);

void update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nt,
                      float   **vx, float    **vy, float    **sxx, float    **syy,
                      float    **sxy, float **pi, float **u, float **uipjp, float **absorb_coeff,
                      float *hc);

void update_s_elastic_VTI_interior ( int * gx, int * gy, int nt,
                                     float **  vx, float **   vy, float **   sxx, float **   syy,float **   sxy,
                                     float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33, GlobVar *gv );

void update_s_elastic_TTI_interior ( int * gx, int * gy, int nt,
                         float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                                    float **   sxx, float **   syy,
                        float **   sxy, float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33,
                                    float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp, GlobVar *gv );

void update_s_elastic_TTI_interior2 ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                        float **  vx, float **   vy, float **   sxx, float **   syy,
                        float **   sxy, float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33,
                                     float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp, float *hc );

void update_s_visc_VTI_interior ( int *gx, int *gy, int nt,
                              float **vx, float **vy, float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                             float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                             float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                             float *bip, float *cip, GlobVar *gv);

void update_s_visc_TTI_interior (int *gx, int *gy, int nt,
                                 float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                                 float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                                  float ** pc11u, float **pc33u, float **pc13u, float ** pc15u, float ** pc35u,
                                 float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                                 float *** pc11d, float ***pc33d, float ***pc13d,
                                 float *** pc15d, float *** pc35d,
                                 float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                                 float *bip, float *cip, GlobVar *gv);

void update_s_elastic_abs(int *gx, int *gy, int nt,
                          float   **vx, float    **vy, float    **sxx, float    **syy,
                          float    **sxy, float **pi, float **u, float **uipjp,
                          float **absorb_coeff, GlobVar *gv);

void update_s_elastic_vti_abs ( int * gx, int * gy,
                               float **  vx, float **   vy, float **   sxx, float **   syy,
                               float **   sxy, float ** pc11, float ** pc55ipjp,
                               float ** pc13, float ** pc33, float ** absorb_coeff, GlobVar *gv );


void update_s_elastic_tti_abs ( int *gx, int *gy, float  **pvxx, float **pvyy, float **pvyx, float **pvxy,
				float **sxx, float **syy, float **sxy, float **pc11, float **pc55ipjp, float **pc13, float **pc33,
				float **pc15, float **pc35, float **pc15ipjp, float **pc35ipjp, float **absorb_coeff);

void update_s_visc_vti_abs ( int *gx, int *gy,
                              float **vx, float **vy, float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                             float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                             float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                             float *bip, float *cip, float ** absorb_coeff, GlobVar *gv );

void update_s_visc_tti_abs ( int *gx, int *gy,
                            float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                            float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                            float ** pc11u, float **pc33u, float **pc13u, float ** pc15u, float ** pc35u,
                           float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                           float *** pc11d, float ***pc33d, float ***pc13d,
                           float *** pc15d, float *** pc35d,
                           float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                           float *bip, float *cip, float ** absorb_coeff, GlobVar *gv );



void update_s_elastic_VTI_PML ( int nx2, int ny2, int * gx, int * gy, 
                            float **  vx, float **   vy, float **   sxx, float **   syy,
                            float **   sxy, float ** pc11, float ** pc13, float **pc33, float ** pc55ipjp,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                            float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, GlobVar *gv );

void update_s_elastic_TTI_PML ( int nx2, int ny2, int * gx, int * gy, 
                               float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                               float **   sxx, float **   syy,
                            float **   sxy,
                               float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33,
                                           float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                               float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, GlobVar *gv );

void update_s_elastic_TTI_PML2 ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                            float **  vx, float **   vy, float **   sxx, float **   syy,
                            float **   sxy,
                               float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33,
                                           float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp, float *hc,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                                float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx );


void update_s_visc_VTI_PML ( int nx2, int ny2, int *gx, int *gy,
                              float **vx, float **vy, float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                             float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                             float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                             float *bip, float *cip,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                            float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, GlobVar *gv );

void update_s_visc_TTI_PML ( int nx2, int ny2, int *gx, int *gy,
                            float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                            float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                            float ** pc11u, float **pc33u, float **pc13u, float ** pc15u, float ** pc35u,
                           float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                           float *** pc11d, float ***pc33d, float ***pc13d,
                           float *** pc15d, float *** pc35d,
                           float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                             float *bip, float *cip,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                            float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, GlobVar *gv );

void update_s_elastic_abs_4(int *gx, int *gy, int nt,
                            float   **vx, float    **vy, float    **sxx, float    **syy,
                            float    **sxy, float **pi, float **u, float **uipjp,
                            float **absorb_coeff, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,
			    float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,
			    float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,
			    float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4, GlobVar *gv);

void update_s_elastic_interior(int *gx, int *gy, int nt,
                               float   **vx, float    **vy, float    **sxx, float    **syy,
                               float    **sxy, float **pi, float **u, float **uipjp,
                               float *hc, GlobVar *gv);

void update_s_elastic_interior_4(int *gx, int *gy, int nt,
                                 float   **vx, float    **vy, float    **sxx, float    **syy,
                                 float    **sxy, float **pi, float **u, float **uipjp, float *hc,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4, GlobVar *gv);

void update_s_elastic_PML(int nx2, int ny2, int *gx, int *gy, int nt,
                          float   **vx, float    **vy, float    **sxx, float    **syy,
                          float    **sxy, float **pi, float **u, float **uipjp,
                          float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                          float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                          float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx, GlobVar *gv);



void update_s_elastic_PML_4(int nx2, int ny2, int *gx, int *gy, int nt,
                            float   **vx, float    **vy, float    **sxx, float    **syy,
                            float    **sxy, float **pi, float **u, float **uipjp,
                            float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                            float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                            float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4, GlobVar *gv);

void update_s_visc_interior(int *gx, int *gy, int nt,
                            float **vx, float **vy, float **sxx, float **syy,
                            float **sxy, float ***r, float *** p, float ***q,
                            float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                            float *cjm, float ***d, float ***e, float ***dip, GlobVar *gv);

void update_s_visc_PML(int nx2, int ny2, int *gx, int *gy,int nt,
                       float   **vx, float    **vy, float    **sxx, float    **syy,
                       float    **sxy, float ***r, float ***p, float ***q, float **fipjp, float **f, float **g,
                       float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip,
                       float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                       float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                       float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx, GlobVar *gv);

void update_s_visc_interior_4(int *gx, int *gy, int nt,
                              float **vx, float **vy, float **sxx, float **syy,
                              float **sxy, float ***r, float *** p, float ***q,
                              float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                              float *cjm, float ***d, float ***e, float ***dip,
                              float *hc,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv);

void update_s_visc_PML_4(int nx2, int ny2, int *gx, int *gy,int nt,
                         float   **vx, float    **vy, float    **sxx, float    **syy,
                         float    **sxy, float ***r, float ***p, float ***q, float **fipjp, float **f, float **g,
                         float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip,
                         float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                         float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                         float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv);

void PML_pro(float *d_x, float *K_x, float *alpha_prime_x, float *a_x, float *b_x,
             float *d_x_half, float *K_x_half, float *alpha_prime_x_half, float *a_x_half, float *b_x_half,
             float *d_y, float *K_y, float *alpha_prime_y, float *a_y, float *b_y,
             float *d_y_half, float *K_y_half, float *alpha_prime_y_half, float *a_y_half, float *b_y_half, GlobVar *gv);


void update_v(int nx1, int nx2, int ny1, int ny2, int nt,
              float   **pvx, float **pvy, float **psxx, float **psyy,
              float **psxy, float **prho, float  **prip, float **prjp,
              float   **srcpos_loc, float **signals, int nsrc, float **absorb_coeff,
              float *hc);

void update_v_abs(int *gx, int *gy, float   **vx, float **vy,
                  float **sxx, float **syy, float **sxy,  float  **rip, float **rjp,
                  float **absorb_coeff, GlobVar *gv);

void update_v_abs_4(int *gx, int *gy, int nt,
                    float **vx, float **vy, float **sxx, float **syy, float **sxy,
                    float **rip, float **rjp, float **absorb_coeff, 
		    float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4, GlobVar *gv);

void update_v_interior(int *gx, int *gy, int nt,
                       float   **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float  **rip, float **rjp,
                       float   **srcpos_loc, float **signals, int nsrc,
                       float *hc, GlobVar *gv);

void update_v_interior_4(int *gx, int *gy, int nt,
                         float   **vx, float **vy, float **sxx, float **syy,
                         float **sxy, float  **rip, float **rjp,
                         float   **srcpos_loc, float **signals, int nsrc,float *hc,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4, GlobVar *gv);

void update_v_PML(int nx2, int ny2,int *gx, int *gy, int nt, float   **vx, float **vy,
                  float **sxx, float **syy, float **sxy, float  **rip, float **rjp,
                  float *K_x, float *a_x, float *b_x, float *K_x_half,
                  float *a_x_half, float *b_x_half, float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half, float **psi_sxx_x, float **psi_syy_y,
                  float **psi_sxy_y, float **psi_syx_x, GlobVar *gv);
void update_v_PML_4(int nx2, int ny2, int *gx, int *gy, int nt, float   **vx, float **vy,
                    float **sxx, float **syy, float **sxy,  float  **rip, float **rjp,
                    float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half,
                    float *b_x_half, float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                    float **psi_sxx_x, float **psi_syy_y, float **psi_sxy_y, float **psi_sxy_x,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4, GlobVar *gv);

void v_derivatives(float **vx, float **vy, float **pvxx, float **pvyy, float **pvyx, float **pvxy, GlobVar *gv);

void wavefield_update_s_el(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                           float **sxx, float **syy, float **pi, float **u, float **uipjp, GlobVar *gv);

void wavefield_update_s_el_vti ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy, float **sxx, float ** syy, float ** pc11, float ** pc55ipjp,
                                float ** pc13, float ** pc33);

void wavefield_update_s_el_tti ( int i, int j,float   **vxx, float  **vyx,float **vxy,float  **vyy,
                                float **sxy, float **sxx, float ** syy, float ** pc11, float ** pc55ipjp,
                                float ** pc13, float ** pc33,
                                float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp);

void wavefield_update_s_el_tti_pml ( int i, int j,float   **vxx, float  **vyx,float **vxy,float  **vyy,
                                float **sxy, float **sxx, float ** syy, float ** pc11, float ** pc55ipjp,
                                float ** pc13, float ** pc33,
                                float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp);



void wavefield_update_s_el_tti2 ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy, float **sxx, float ** syy, float ** pc11,                                 float ** pc55ipjp,
                            float ** pc13, float ** pc33,
                                 float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp);

void wavefield_update_s_visc(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                             float **sxx, float **syy, float ***r, float ***p,
                             float ***q,float **fipjp, float **f, float **g, float *bip,
                             float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip, GlobVar *gv);

void wavefield_update_s_visc_VTI ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy,
                              float **sxy, float **sxx, float ** syy, float ***pp, float ***pr,float ***pq,
                              float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                              float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                              float *bip, float *cip, GlobVar *gv);

void wavefield_update_s_visc_TTI ( int i, int j,float  **vxx, float **vyx,float **vxy,float **vyy,
                              float **sxy, float **sxx, float ** syy, float ***p, float ***r,float ***q,
                                  float ** pc11u, float **pc33u, float **pc13u, float ** pc15u, float ** pc35u,
                                 float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                                 float *** pc11d, float ***pc33d, float ***pc13d,
                                 float *** pc15d, float *** pc35d,
                                 float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                                  float *bip, float *cip, GlobVar *gv);

void wavefield_update_s_visc_4(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                               float **sxx, float **syy, float ***r, float ***p,
                               float ***q,float **fipjp, float **f, float **g, float *bip,
                               float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv);

void wavefield_update_v(int i, int j,float   sxx_x, float  sxy_x,float sxy_y,float  syy_y, float **vx,
                        float **vy, float **rip, float **rjp, GlobVar *gv);
void wavefield_update_v_4(int i, int j,float   sxx_x, float  sxy_x,float sxy_y,float  syy_y, float **vx,
                          float **vy, float **rip, float **rjp,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4, GlobVar *gv);

void wavefield_update_s_el_4(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                             float **sxx, float **syy, float **pi, float **u, float **uipjp,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4, GlobVar *gv);

float **wavelet(float **srcpos_loc, int nsrc, GlobVar *gv);

void writebufs(float **sxx, float **syy,
               float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
               float **buffertop_to_bot, float **bufferbot_to_top);

void writebufv(float **vx, float **vy,
               float **bufferlef_to_rig, float **bufferrig_to_lef,
               float **buffertop_to_bot, float **bufferbot_to_top);

void write_par(GlobVar *gv);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char *modfile, float **array, int format, GlobVar *gv);

void zero_elastic(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx, float **syy, float **sxy);

void zero_elastic_4(int nx1, int nx2, int ny1, int ny2, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void zero_visco_4(int nx1, int nx2, int ny1, int ny2, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv);

void zero_visc(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx, float **syy, float **sxy, float *** pr, float *** pp, float *** pq, GlobVar *gv);

void zero_PML_elastic(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx,
                      float **syy, float **sxy,
                      float **psi_sxx_x, float **psi_sxy_x, float **psi_vxx, float **psi_vyx, float **psi_syy_y, float **psi_sxy_y, float **psi_vyy, float **psi_vxy,float **psi_vxxs, GlobVar *gv);

void zero_PML_visc(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx,
                   float **syy, float **sxy,
                   float **psi_sxx_x, float **psi_sxy_x, float **psi_vxx, float **psi_vyx, float **psi_syy_y, float **psi_sxy_y, float **psi_vyy, float **psi_vxy, float **psi_vxxs,
                   float ***pr, float ***pp, float ***pq, GlobVar *gv);

/* declaration of functions for parser*/

/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(const char* input_file,char **varname_list,char **value_list);

void print_objectlist_screen(int number_readobject,char **varname_list,char **value_list);

int count_occure_charinstring(char *stringline, char teststring[]);

void copy_str2str_uptochar(char *string_in, char *string_out, char teststring[]);

int get_int_from_objectlist(char *string_in, int number_readobject, int *int_buffer,
                            char **varname_list,char **value_list, int *used_list);

int get_float_from_objectlist(char *string_in, int number_readobject, float *double_buffer,
                              char **varname_list,char **value_list, int *used_list);

int get_string_from_objectlist(char *string_in, int number_readobject, char *string_buffer,
                               char **varname_list,char **value_list, int *used_list);

int is_string_blankspace(char *string_in);

void remove_blankspaces_around_string(char *string_in);

void add_object_tolist(char *string_name,char *string_value, int *number_read_object,
                       char **varname_list,char **value_list);

/* utility functions */
void dt_mult(int nx, int ny, float dt, float  **  a );
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

#endif
