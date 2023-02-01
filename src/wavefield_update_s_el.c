
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

/* void wavefield_update_s_el(int i, int j, float vxx, float vyx, float vxy, float vyy, float **sxy,
                           float **sxx, float **syy, float **pi, float **u, float **uipjp, GlobVar *gv) */

void wavefield_update_s_el(int i, int j, float vxx, float vyx, float vxy, float vyy, MemModel * mpm, MemWavefield * mpw,
                           GlobVar * gv)
{
    float fipjp = mpm->puipjp[j][i] * gv->DT;
    float f = mpm->pu[j][i] * gv->DT;
    float g = mpm->ppi[j][i] * gv->DT;

    /* Update  */
    mpw->psxy[j][i] += fipjp * (vyx + vxy);
    mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
    mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
}
