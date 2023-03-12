
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
 * memm_struct.h - variables used for memory initiallisation.
 * ----------------------------------------------------------------------*/

#ifndef MEMM_STRUCT_H_INCLUDED
#define MEMM_STRUCT_H_INCLUDED

typedef struct {

    float **prho;
    float **prip;
    float **prjp;
    float **ppi;

    float **pu;
    float **puipjp;

    float **absorb_coeff;

    float **ptaus;
    float **ptaup;
    float *peta;
    float **ptausipjp;
    float **fipjp;

    float ***dip;
    float *bip;
    float *bjm;
    float *cip;
    float *cjm;
    float ***d;
    float ***e;
    float **f;
    float **g;

    float **pc11;
    float **pc33;
    float **pc13;
    float **pc55;
    float **pc15;
    float **pc35;
    float **pc55ipjp;
    float **pc15ipjp;
    float **pc35ipjp;
    float **pc11u;
    float **pc33u;
    float **pc13u;
    float **pc55u;
    float **pc15u;
    float **pc35u;
    float **pc55ipjpu;
    float **pc15ipjpu;
    float **pc35ipjpu;
    float ***pc11d;
    float ***pc33d;
    float ***pc13d;
    float ***pc55d;
    float ***pc15d;
    float ***pc35d;
    float ***pc55ipjpd;
    float ***pc15ipjpd;
    float ***pc35ipjpd;

    float **ptau11;
    float **ptau33;
    float **ptau13;
    float **ptau55;
    float **ptau15;
    float **ptau35;
    float **ptau55ipjp;
    float **ptau15ipjp;
    float **ptau35ipjp;

    /* PML variables */
    float *d_x;
    float *K_x;
    float *alpha_prime_x;
    float *a_x;
    float *b_x;
    float *d_x_half;
    float *K_x_half;
    float *alpha_prime_x_half;
    float *a_x_half;
    float *b_x_half;
    float *d_y;
    float *K_y;
    float *alpha_prime_y;
    float *a_y;
    float *b_y;
    float *d_y_half;
    float *K_y_half;
    float *alpha_prime_y_half;
    float *a_y_half;
    float *b_y_half;

/*    float **psi_sxx_x;
    float **psi_syy_y;
    float **psi_sxy_y;
    float **psi_sxy_x;
    float **psi_vxx;
    float **psi_vyy;
    float **psi_vxy;
    float **psi_vyx;
    float **psi_vxxs; */

} MemModel;

#endif
