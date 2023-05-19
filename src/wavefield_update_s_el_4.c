
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

void wavefield_update_s_el_4(int i, int j, MemModel * mpm, MemWavefield * mpw, GlobVar *gv)
{
    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;
    float sumxx = 0.0f, sumyy = 0.0f, sumxy = 0.0f, sumyx = 0.0f;
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;

    float dhi = 1.0 / gv->DH;
    /*float fipjp = mpm->puipjp[j][i] * gv->DT;
    float f = mpm->pu[j][i] * gv->DT;
    float g = mpm->ppi[j][i] * gv->DT;*/

    // Save 
    mpw->vxx_1[j][i] = mpw->pvxx[j][i] * gv->DH;
    mpw->vyy_1[j][i] = mpw->pvyy[j][i] * gv->DH;
    mpw->vxy_1[j][i] = mpw->pvxy[j][i] * gv->DH;
    mpw->vyx_1[j][i] = mpw->pvyx[j][i] * gv->DH;

    // Calculate Adams-Bashforth stuff
    sumxx = c1 * mpw->vxx_1[j][i] + c2 * mpw->vxx_2[j][i] + c3 * mpw->vxx_3[j][i] + c4 * mpw->vxx_4[j][i];
    sumyy = c1 * mpw->vyy_1[j][i] + c2 * mpw->vyy_2[j][i] + c3 * mpw->vyy_3[j][i] + c4 * mpw->vyy_4[j][i];
    sumxy = c1 * mpw->vxy_1[j][i] + c2 * mpw->vxy_2[j][i] + c3 * mpw->vxy_3[j][i] + c4 * mpw->vxy_4[j][i];
    sumyx = c1 * mpw->vyx_1[j][i] + c2 * mpw->vyx_2[j][i] + c3 * mpw->vyx_3[j][i] + c4 * mpw->vyx_4[j][i];

    /* The stress updates are stored in internal variables due to extensions necessary for FWI */
    /* calculate stress component update */
    u1 = mpm->fipjp[j][i] * (sumxy + sumyx) * dhi;
    u2 = ((mpm->g[j][i] * (sumxx + sumyy)) - (2.0 * mpm->f[j][i] * sumyy)) * dhi;
    u3 = ((mpm->g[j][i] * (sumxx + sumyy)) - (2.0 * mpm->f[j][i] * sumxx)) * dhi;
    
    /* updating components of the stress tensor */
    mpw->psxy[j][i] += u1;
    mpw->psxx[j][i] += u2;
    mpw->psyy[j][i] += u3;
}
