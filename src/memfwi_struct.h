
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

/*------------------------------------------------------------------------
 * memfwi_struct.h - variables used for memory initiallisation of inversion parameters.
 *
 * ----------------------------------------------------------------------*/

#ifndef MEMINV_STRUCT_H_INCLUDED
#define MEMINV_STRUCT_H_INCLUDED

typedef struct {

    float *L2_hist;                     // vector for abort criterion
    
    int *DTINV_help;                    // time sample increment used for inversion

    float **ux;
    float **uy;
    float **uxy;
    float **uyx;
    //float **u;
    float **uttx;
    float **utty;

    float **pvxp1;
    float **pvyp1;
    float **pvxm1;
    float **pvym1;
    
    /* Model allocations */
    float **Vp0;
    float **Vs0;
    float **Rho0;

    float **vpmat;
    
    float **prhonp1;
    float **pripnp1;
    float **prjpnp1;
    float **ppinp1;
    float **punp1;

    float ***forward_prop_x;
    float ***forward_prop_y;
    float **gradg;
    float **gradp;
    float ***forward_prop_rho_x;
    float ***forward_prop_rho_y;
    float **gradg_rho;
    float **gradp_rho;
    float ***forward_prop_u;
    float **gradg_u;
    float **gradp_u;
    //float ***forward_prop_p;

    float **waveconv;
    float **waveconv_lam;
    float **waveconv_shot;
    float **waveconvtmp;
    float **wcpart;
    float **wavejac;
    float **waveconv_rho;
    float **waveconv_rho_s;
    float **waveconv_rho_shot;
    float **waveconv_u;
    float **waveconv_mu;
    float **waveconv_u_shot;
    
    /* for Wolfe condition*/

    float **waveconv_old;
    float **waveconv_u_old;
    float **waveconv_rho_old;
    float **waveconv_up;
    float **waveconv_u_up;
    float **waveconv_rho_up;
    /*int steplength_search = 0;
    int FWI_run = 1;
    int gradient_optimization = 1;
    float alpha_SL_min = 0;
    float alpha_SL_max = 0;
    float alpha_SL = 1.0;
    float alpha_SL_old;
    float L2_SL_old = 0;
    float L2_SL_new = 0;
    float c1_SL = 1e-4;
    float c2_SL = 0.9;
    int wolfe_status;
    int wolfe_sum_FWI = 0;
    int wolfe_found_lower_L2 = 0;
    float alpha_SL_FS;
    float L2_SL_FS;
    int use_wolfe_failsafe = 0;
    int wolfe_SLS_failed = 0; */

    /* Variables for L-BFGS */
    float **s_LBFGS;            //
    float **y_LBFGS;            //
    float *rho_LBFGS;           //

    
    /* Variables for energy weighted gradient */
    float **Ws;                 // total energy of the source wavefield
    float **Wr;                 // total energy of the receiver wavefield
    float **We;                 // total energy of source agv->ND receiver wavefield
    float **We_sum;
    //float We_sum_max1;
    //float We_max;

    float **fulldata;
    float **fulldata_vx;
    float **fulldata_vy;
    float **fulldata_vz;
    float **fulldata_p;
    float **fulldata_curl;
    float **fulldata_div;
    
    float **sectionpnp1;
    float **sectionpn;
    float **sectionpdata;
    float **sectionpdiff;
    float **sectionpdiffold;
    float **sectionvxdata;
    float **sectionvxdiff;
    float **sectionvxdiffold;
    float **sectionvydata;
    float **sectionvydiff;
    float **sectionvydiffold;

    float **sectionread;
    float **sectionvy_conv;
    float **sectionvy_obs;
    float **sectionvx_conv;
    float **sectionvx_obs;
    float **sectionp_conv;
    float **sectionp_obs;
    float *source_time_function;

    float **taper_coeff;
    float *L2t;
    float *epst1;
    float *epst2;
    float *epst3;
    float *picked_times;

} MemInv;

#endif
