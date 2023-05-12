
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

void wavefield_update_s_el_vti(int i, int j, MemModel * mpm, MemWavefield * mpw, MemInv *minv, GlobVar * gv)
{
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;
    
    /* calculate stress component update */
    u1 = (mpm->pc55ipjp[j][i] * (mpw->pvxy[j][i] + mpw->pvyx[j][i]));
    u2 = ((mpm->pc11[j][i] * mpw->pvxx[j][i]) + (mpm->pc13[j][i] * mpw->pvyy[j][i]));
    u3 = ((mpm->pc13[j][i] * mpw->pvxx[j][i]) + (mpm->pc33[j][i] * mpw->pvyy[j][i]));

    /* updating components of the stress tensor */
    mpw->psxy[j][i] += u1;
    mpw->psxx[j][i] += u2;
    mpw->psyy[j][i] += u3;
    
    /* updating components of the gradient */
    if (gv->MODE == FWI) {
        minv->uxy[j][i] = u1 / gv->DT;
        minv->ux[j][i] = u2 / gv->DT;
        minv->uy[j][i] = u3 / gv->DT;
    }
}
