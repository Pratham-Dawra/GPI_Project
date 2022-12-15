
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

/*
 * Update Function of the stress-Wavefields in the elastic case
 */

#include "fd.h"

void wavefield_update_s_el_tti(int i, int j, float **vxx, float **vyx, float **vxy, float **vyy, float **sxy,
                               float **sxx, float **syy, float **pc11, float **pc55ipjp, float **pc13, float **pc33,
                               float **pc15, float **pc35, float **pc15ipjp, float **pc35ipjp)
{
    float vxxipjp = 0.25 * (vxx[j][i] + vxx[j + 1][i] + vxx[j][i + 1] + vxx[j + 1][i + 1]);
    float vyyipjp = 0.25 * (vyy[j][i] + vyy[j + 1][i] + vyy[j][i + 1] + vyy[j + 1][i + 1]);   
    float vij = (0.25 * (vyx[j][i] + vyx[j - 1][i] + vyx[j][i - 1] + vyx[j - 1][i - 1])) + 
                (0.25 * (vxy[j][i] + vxy[j - 1][i] + vxy[j][i - 1] + vxy[j - 1][i - 1]));
    float v = vxy[j][i] + vyx[j][i];

    /* Update  */
    sxx[j][i] += ((pc11[j][i] * vxx[j][i]) + (pc13[j][i] * vyy[j][i]) + (pc15[j][i] * vij));
    syy[j][i] += ((pc13[j][i] * vxx[j][i]) + (pc33[j][i] * vyy[j][i]) + (pc35[j][i] * vij));

    sxy[j][i] += ((pc55ipjp[j][i] * v) + (pc15ipjp[j][i] * vxxipjp) + (pc35ipjp[j][i] * vyyipjp));
}
