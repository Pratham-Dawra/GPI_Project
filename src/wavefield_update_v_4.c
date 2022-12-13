
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
 * Update Function of the particle velocity wavefields 
 */

#include "fd.h"

void wavefield_update_v_4(int i, int j, float sxx_x, float sxy_x, float sxy_y, float syy_y, float **vx,
                          float **vy, float **rip, float **rjp, float **svx_1, float **svx_2, float **svx_3,
                          float **svx_4, float **svy_1, float **svy_2, float **svy_3, float **svy_4, GlobVar * gv)
{
    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;

    float dtdh = gv->DT / gv->DH;

    // Save derivations
    svx_1[j][i] = sxx_x + sxy_y;
    svy_1[j][i] = sxy_x + syy_y;

    vx[j][i] += (c1 * svx_1[j][i] + c2 * svx_2[j][i] + c3 * svx_3[j][i] + c4 * svx_4[j][i]) * dtdh * rip[j][i];
    vy[j][i] += (c1 * svy_1[j][i] + c2 * svy_2[j][i] + c3 * svy_3[j][i] + c4 * svy_4[j][i]) * dtdh * rjp[j][i];
}
