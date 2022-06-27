/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/* $Id: wavefield_update_s_el.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the elastic case*/

#include "fd.h"

void wavefield_update_s_el_tti ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy, float **sxx, float ** syy, float ** pc11, float ** pc55ipjp,
                            float ** pc13, float ** pc33,
                            float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp)
{

    float c55ipjp, c33, c11, c13;
    float c15ipjp, c35ipjp, c15, c35, v;


    c55ipjp=pc55ipjp[j][i];
    c35ipjp=pc35ipjp[j][i];
    c15ipjp=pc15ipjp[j][i];

    c11=pc11[j][i];
    c33=pc33[j][i];
    c13=pc13[j][i];
    c15=pc15[j][i];
    c35=pc35[j][i];

    v=vxy+vyx;

	/*Update  */
    sxy[j][i]+= ( (c55ipjp*v) + (c15ipjp*vxx)+(c35ipjp*vyy));
    sxx[j][i]+= ( (c11* vxx)+ (c13*vyy) + (c15*v));
    syy[j][i]+= ( (c13* vxx)+ (c33*vyy) + (c35*v));
    

}
