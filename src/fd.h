
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

/*----------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD program sofi2D
 * ---------------------------------------------------------------------*/

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
#include "globvar_struct_fwi.h"
#include "acq_struct.h"
#include "memm_struct.h"
#include "memfwi_struct.h"
#include "memw_struct.h"
#include "perform_struct.h"
#include "util.h"
#include "macros.h"

/* declaration of functions */
void abs_update_s(int i, int j, MemModel *mpm, MemWavefield *mpw);

// void abs_update_s_ac1(int i, int j, MemModel *mpm, MemWavefield *mpw);

void abs_update_s_ac2(int i, int j, MemModel *mpm, MemWavefield *mpw);

void abs_update_v(int i, int j, MemModel *mpm, MemWavefield *mpw);

void absorb(float **absorb_coeff, GlobVar *gv);

int acq_read(AcqVar *acq, GlobVar *gv);

void apply_workflow(int iter, GlobVar *gv, GlobVarInv *vinv);

void av_mat(float **pi, float **u, float **ppijm, float **puip, float **pujm);

void av_mue(float **u, float **uipjp, GlobVar *gv);

void av_rho(float **rho, float **rip, float **rjp, GlobVar *gv);

void av_tau(float **taus, float **tausipjp, GlobVar *gv);

void calc_conv(int hin, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void calc_envelope(float **datatrace, float **envelope, int ntr, int ns);

void calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold,
              int sws, int swstestshot, int ishot, int iter, AcqVar *acq, GlobVar *gv, GlobVarInv *vinv);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);

void check_fs(GlobVar *gv, GlobVarInv *vinv);

void checkfd(float *hc, float **srcpos, int nsrc, int **recpos, GlobVar *gv);

void cpml_update_s_x(int i, int j, int h1, int h2, MemModel *mpm, MemWavefield *mpw);

void cpml_update_s_x_ac(int i, int j, int h1, int h2, MemModel *mpm, MemWavefield *mpw);

void cpml_update_s_y(int i, int j, int h1, int h2, MemModel *mpm, MemWavefield *mpw);

void cpml_update_s_y_ac(int i, int j, int h1, int h2, MemModel *mpm, MemWavefield *mpw);

void cpml_update_v_x(int i, int j, float *sxx_x, float *sxy_x, MemModel *mpm, MemWavefield *mpw);

void cpml_update_v_x_ac(int i, int j, float *sxx_x, MemModel * mpm, MemWavefield * mpw);

void cpml_update_v_y(int i, int j, float *sxy_y, float *syy_y, MemModel *mpm, MemWavefield *mpw);

void cpml_update_v_y_ac(int i, int j, float *syy_y, MemModel * mpm, MemWavefield * mpw);

void eprecond(MemWavefield *mpw, float **W, GlobVar *gv);

void eqsource(int nt, AcqVar *acq, MemWavefield *mpw, GlobVar *gv);

void exchange_v(float **vx, float **vy, MemWavefield *mpw, GlobVar *gv);

void exchange_s(MemWavefield *mpw, GlobVar *gv);

void exchange_par(GlobVar *gv, GlobVarInv *vinv);

void filter_frequencies(GlobVarInv *vinv);

void freemem(MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void freemem_model(MemModel *mpm, GlobVar *gv);

void freemem_wavefield(MemWavefield *mpw, GlobVar *gv);

int get_fd_order(GlobVar *gv);

const char *get_weq_verbose(WEQTYPE wt);

const char *get_mode_verbose(RUNMODE md);

float *holbergcoeff(GlobVar *gv);

void initfd(GlobVar *gv);

void init_grad(int iter, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void init_inv(int iter, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void initmem(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void initmem_fwi(MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void initmem_model(MemModel *mpm, GlobVar *gv);

void initmem_wavefield(MemWavefield *mpw, GlobVar *gv);

void initproc(GlobVar *gv);

int initsrc(int ishot, int nshots, AcqVar *acq, GlobVar *gv);

void inseis(int ishot, int sws, float **section, GlobVar *gv, GlobVarInv *vinv);

void inversion(int iter, int ishot, int snapcheck, float *hc, AcqVar *acq, MemModel *mpm, MemWavefield *mpw,
               MemInv *minv, GlobVar *gv, GlobVarInv *vinv, Perform *perf);

void lbfgs_reset(int iter, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void lbfgs_core(int iter, float *q_LBFGS, float *alpha_LBFGS, float *r_LBFGS, MemInv *minv, GlobVar *gv,
                GlobVarInv *vinv);

void lbfgs(float **grad_vs, float **grad_rho, float **grad_vp, float Vs_avg, float rho_avg, float Vp_avg, int iter,
           MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void model_elastic(MemModel *mpm, GlobVar *gv);

void model_elastic_VTI(MemModel *mpm, GlobVar *gv);

void model_elastic_TTI(MemModel *mpm, GlobVar *gv);

void model_visco(MemModel *mpm, GlobVar *gv);

void model_visco_vti(MemModel *mpm, GlobVar *gv);

void model_visco_tti(MemModel *mpm, GlobVar *gv);

void matcopy(float **prho, float **ppi, float **pu, float **ptaup, float **ptaus, GlobVar *gv);

void matcopy_elastic(float **prho, float **ppi, float **pu, GlobVar *gv);

void matcopy_ani(float **rho, float **c11, float **c15, float **c13,
                 float **c35, float **c33, float **c55, float **taus, float **taup);

void merge(int nsnap, int type, int SNAPIDX[][5], GlobVar *gv);

void mergemod(const char *modfile, int format, GlobVar *gv);

void outseis_glob(FILE *fpdata, float **section,
                  int **recpos, int ntr, float **srcpos_loc, int ns, int seis_form, int ishot, int comp, GlobVar *gv);

void output_gradient(int ishot, int nshots, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void par_mult_dt(float **pi, float **u, float **uipjp);

void PML_pro(MemModel *mpm, GlobVar *gv);

void prepare_update_s_ac(MemModel *mpm, GlobVar *gv);

void prepare_update_s_el(MemModel *mpm, GlobVar *gv);

void prepare_update_s_el_4(MemModel *mpm, GlobVar *gv);

void prepare_update_s_visc(MemModel *mpm, GlobVar *gv);

void prepare_update_s_visc_4(MemModel *mpm, GlobVar *gv);

void prepare_update_s_vti(MemModel *mpm, GlobVar *gv);

void prepare_update_s_tti(MemModel *mpm, GlobVar *gv);

void prepmod(MemModel *mpm, GlobVar *gv);

void psource(int nt, float **srcpos, float **signal, int nsrc, MemWavefield *mpw, GlobVar *gv);

void psource_rsg(int nt, float **sxx, float **syy, float **srcpos_loc, float **signals, int nsrc, GlobVar *gv);

float readdsk(FILE *fp_in, int format);

void readbufv(float **vx, float **vy,
              float **bufferlef_to_rig, float **bufferrig_to_lef, float **buffertop_to_bot, float **bufferbot_to_top);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float **vx, float **vy, float **sxx, float **syy, float **sxy, GlobVar *gv);

void read_par_json(const char *fileinp, GlobVar *gv, GlobVarInv *vinv);

int read_par_json_fwi(int number_readobjects, char **varname_list, char **value_list, int *used_list, int fserr,
                      GlobVar *gv, GlobVarInv *vinv);

void read_workflow(GlobVar *gv, GlobVarInv *vinv);

void readmod(MemModel *mpm, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void readmod_acoustic(MemModel *mpm, GlobVar *gv);

void readmod_acoustic_vti(MemModel *mpm, GlobVar *gv);

void readmod_acoustic_tti(MemModel *mpm, GlobVar *gv);

void readmod_elastic(MemModel *mpm, MemInv *minv, GlobVar *gv);

void readmod_elastic_vti(MemModel *mpm, MemInv *minv, GlobVar *gv);

void readmod_elastic_tti(MemModel *mpm, MemInv *minv, GlobVar *gv);

void readmod_visco(MemModel *mpm, MemInv *minv, GlobVar *gv);

void readmod_visco_vti(MemModel *mpm, MemInv *minv, GlobVar *gv);

void readmod_visco_tti(MemModel *mpm, MemInv *minv, GlobVar *gv);

int **receiver(GlobVar *gv, int *topo);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float **vx, float **vy, float **sxx, float **syy, float **sxy, GlobVar *gv);

void saveseis(int ishot, AcqVar *acq, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void saveseis_fwi(int ishot, AcqVar *acq, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void saveseis_glob(float **sectiondata, int **recpos, float **srcpos, int ishot, int ns, int sectiondatatype,
                   GlobVar *gv);

int *scan_topo(GlobVar *gv);

void set_fd_order(int new_order, GlobVar *gv);

void seismo_ssg(int lsamp, int **recpos, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void snap(int nt, int nsnap, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void snap_store(int nt, int hin, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void snapmerge(int nsnap);

void sources(AcqVar *acq, GlobVar *gv, int *topo);

int **splitrec(int **recpos, int *recswitch, GlobVar *gv);

float **splitsrc(float **srcpos, int *nsrc_loc, int nsrc, GlobVar *gv);

void subgrid_bounds(int nx1, int nx2, int ny1, int ny2, GlobVar *gv);

void surface(int ndepth, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void surface_elastic(int ndepth, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void timedomain_filt(float **data, int method, int ntr, GlobVar *gv, GlobVarInv *vinv);

void time_loop(int iter, int ishot, int snapcheck, float *hc, float **srcpos_loc, float **signalx,
               float **signaly, int nsrc, int sw, AcqVar *acq, MemModel *mpm, MemWavefield *mpw,
               MemInv *minv, GlobVar *gv, GlobVarInv *vinv, Perform *perf);

/* void update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nt,
                      float **vx, float **vy, float **sxx, float **syy,
                      float **sxy, float **pi, float **u, float **uipjp, float **absorb_coeff, float *hc); */

void update_s_acoustic_abs(int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_s_acoustic_interior(int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_s_elastic_abs(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_abs_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_interior_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_acoustic_PML(int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_s_acoustic_vti_interior(int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_s_elastic_PML(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_PML_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_vti_abs(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_vti_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_vti_PML(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_tti_abs(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_elastic_tti_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

/* void update_s_elastic_TTI_interior2(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                                    float **vx, float **vy, float **sxx, float **syy,
                                    float **sxy, float **pc11, float **pc55ipjp, float **pc13, float **pc33,
                                    float **pc15, float **pc35, float **pc15ipjp, float **pc35ipjp, float *hc); */

void update_s_elastic_tti_PML(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

/* void update_s_elastic_TTI_PML2(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                               float **vx, float **vy, float **sxx, float **syy,
                               float **sxy,
                               float **pc11, float **pc55ipjp, float **pc13, float **pc33,
                               float **pc15, float **pc35, float **pc15ipjp, float **pc35ipjp, float *hc,
                               float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                               float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                               float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx); */

void update_s_visc_abs(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_abs_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_interior_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_PML(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_PML_4(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_vti_abs(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_vti_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_vti_PML(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_tti_abs(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_tti_interior(int nt, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void update_s_visc_tti_PML(MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

/* void update_v(int nx1, int nx2, int ny1, int ny2, int nt,
              float **pvx, float **pvy, float **psxx, float **psyy,
              float **psxy, float **prho, float **prip, float **prjp,
              float **srcpos_loc, float **signals, int nsrc, float **absorb_coeff, float *hc); */

void update_v_abs(int sw, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void update_v_abs_4(int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_v_interior(int nt, float **srcpos_loc, float **signals, float **signals1, int nsrc, int sw,
                       MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void update_v_interior_4(int nt, float **srcpos_loc, float **signals, int nsrc, float *hc,
                         MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void update_v_PML(int nx2, int ny2, int nt, int sw, MemModel *mpm, MemWavefield *mpw, MemInv *minv,
                  GlobVar *gv, GlobVarInv *vinv);

void update_v_PML_4(int nx2, int ny2, int nt, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void v_derivatives(MemWavefield *mpw, GlobVar *gv);

void wavefield_update_s_ac(int i, int j, MemModel *mpm, MemWavefield *mpw);

void wavefield_update_s_ac_vti(int i, int j, MemModel *mpm, MemWavefield *mpw);

void wavefield_update_s_el(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_el_4(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_el_vti(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_el_tti(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_visc(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_visc_4(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_visc_VTI(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_s_visc_TTI(int i, int j, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv);

void wavefield_update_v(int i, int j, int sw, float sxx_x, float sxy_x, float sxy_y, float syy_y, MemModel *mpm,
                        MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

void wavefield_update_v_ac(int i, int j, float sxx_x, float sxx_y, MemModel *mpm, MemWavefield *mpw, GlobVar *gv);

void wavefield_update_v_4(int i, int j, float sxx_x, float sxy_x, float sxy_y, float syy_y, MemModel *mpm,
                          MemWavefield *mpw, GlobVar *gv);

void wavelet(AcqVar *acq, GlobVar *gv);

void writebufs(float **sxx, float **syy,
               float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
               float **buffertop_to_bot, float **bufferbot_to_top);

void writebufv(float **vx, float **vy,
               float **bufferlef_to_rig, float **bufferrig_to_lef, float **buffertop_to_bot, float **bufferbot_to_top);

void write_par(GlobVar *gv, GlobVarInv *vinv);

void write_par_fwi(GlobVar *gv, GlobVarInv *vinv);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(const char *modfile, float **array, int format, const GlobVar *gv);

/*void zero_elastic(int j, int i, MemWavefield *mpw);

void zero_acoustic(int j, int i, MemWavefield *mpw);

void zero_elastic_4(int j, int i, MemWavefield *mpw);

void zero_visco_4(int j, int i, int l, MemWavefield *mpw);

void zero_visco(int j, int i, int l, MemWavefield *mpw);

void zero_PML_x(int j, int i, MemWavefield *mpw);

void zero_PML_y(int j, int i, MemWavefield *mpw);*/

void zero_wavefield(int iter, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv);

/* declaration of functions for json parser in json_parser.c */

int read_objects_from_intputfile(const char *input_file, char **varname_list, char **value_list);

void print_objectlist_screen(int number_readobject, char **varname_list, char **value_list);

int count_occure_charinstring(char *stringline, char teststring[]);

void copy_str2str_uptochar(char *string_in, char *string_out, char teststring[]);

int get_int_from_objectlist(char *string_in, int number_readobject, int *int_buffer,
                            char **varname_list, char **value_list, int *used_list);

int get_float_from_objectlist(char *string_in, int number_readobject, float *double_buffer,
                              char **varname_list, char **value_list, int *used_list);

int get_string_from_objectlist(char *string_in, int number_readobject, char *string_buffer,
                               char **varname_list, char **value_list, int *used_list);

int is_string_blankspace(char *string_in);

void remove_blankspaces_around_string(char *string_in);

void add_object_tolist(char *string_name, char *string_value, int *number_read_object,
                       char **varname_list, char **value_list);

#endif
