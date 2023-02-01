
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
 * mem_struct.h - variables used for memory initiallisation.
 *
 * ----------------------------------------------------------------------*/

#ifndef MEMW_STRUCT_H_INCLUDED
#define MEMW_STRUCT_H_INCLUDED

typedef struct {

    float **psxx;
    float **psxy;
    float **psyy;
    float **pvx;
    float **pvy;

    float **vxx_1;
    float **vxx_2;
    float **vxx_3;
    float **vxx_4;

    float **vyy_1;
    float **vyy_2;
    float **vyy_3;
    float **vyy_4;

    float **vxy_1;
    float **vxy_2;
    float **vxy_3;
    float **vxy_4;

    float **vyx_1;
    float **vyx_2;
    float **vyx_3;
    float **vyx_4;

    float **svx_1;
    float **svx_2;
    float **svx_3;
    float **svx_4;

    float **svy_1;
    float **svy_2;
    float **svy_3;
    float **svy_4;

    float ***pr;
    float ***pr_2;
    float ***pr_3;
    float ***pr_4;

    float ***pp;
    float ***pp_2;
    float ***pp_3;
    float ***pp_4;

    float ***pq;
    float ***pq_2;
    float ***pq_3;
    float ***pq_4;

    float **pvxx;
    float **pvyy;
    float **pvyx;
    float **pvxy;

    /* buffer arrays in which the wavefield information to be exchanged between neighboring PEs is stored */
    float **bufferlef_to_rig;
    float **bufferrig_to_lef;
    float **buffertop_to_bot;
    float **bufferbot_to_top;

    /* PML variables */
    float **psi_sxx_x;
    float **psi_syy_y;
    float **psi_sxy_y;
    float **psi_sxy_x;
    float **psi_vxx;
    float **psi_vyy;
    float **psi_vxy;
    float **psi_vyx;
    float **psi_vxxs;

} MemWavefield;

#endif
