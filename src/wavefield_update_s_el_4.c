
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

void wavefield_update_s_el_4(int i, int j, float vxx, float vyx, float vxy, float vyy, float **sxy,
                             float **sxx, float **syy, float **pi, float **u, float **uipjp, float **vxx_1,
                             float **vxx_2, float **vxx_3, float **vxx_4, float **vyy_1, float **vyy_2, float **vyy_3,
                             float **vyy_4, float **vxy_1, float **vxy_2, float **vxy_3, float **vxy_4, float **vyx_1,
                             float **vyx_2, float **vyx_3, float **vyx_4, GlobVar *gv)
{
    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;
    float sumxx = 0.0f, sumyy = 0.0f, sumxy = 0.0f, sumyx = 0.0f;

    float dhi = 1.0 / gv->DH;
    float fipjp = uipjp[j][i] * gv->DT;
    float f = u[j][i] * gv->DT;
    float g = pi[j][i] * gv->DT;

    // Save 
    vxx_1[j][i] = vxx * gv->DH;
    vyy_1[j][i] = vyy * gv->DH;
    vxy_1[j][i] = vxy * gv->DH;
    vyx_1[j][i] = vyx * gv->DH;

    // Calculate Adams-Bashforth stuff
    sumxx = c1 * vxx_1[j][i] + c2 * vxx_2[j][i] + c3 * vxx_3[j][i] + c4 * vxx_4[j][i];
    sumyy = c1 * vyy_1[j][i] + c2 * vyy_2[j][i] + c3 * vyy_3[j][i] + c4 * vyy_4[j][i];
    sumxy = c1 * vxy_1[j][i] + c2 * vxy_2[j][i] + c3 * vxy_3[j][i] + c4 * vxy_4[j][i];
    sumyx = c1 * vyx_1[j][i] + c2 * vyx_2[j][i] + c3 * vyx_3[j][i] + c4 * vyx_4[j][i];

    // Update stress
    sxy[j][i] += ((fipjp * (sumxy + sumyx))) * dhi;
    sxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
    syy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
}
