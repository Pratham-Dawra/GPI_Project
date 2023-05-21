
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

void wavefield_update_s_el_tti(int i, int j, MemModel * mpm, MemWavefield * mpw)
{
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;
    float vxxipjp = 0.25 * (mpw->pvxx[j][i] + mpw->pvxx[j + 1][i] + mpw->pvxx[j][i + 1] + mpw->pvxx[j + 1][i + 1]);
    float vyyipjp = 0.25 * (mpw->pvyy[j][i] + mpw->pvyy[j + 1][i] + mpw->pvyy[j][i + 1] + mpw->pvyy[j + 1][i + 1]);
    float vij = (0.25 * (mpw->pvyx[j][i] + mpw->pvyx[j - 1][i] + mpw->pvyx[j][i - 1] + mpw->pvyx[j - 1][i - 1])) +
        (0.25 * (mpw->pvxy[j][i] + mpw->pvxy[j - 1][i] + mpw->pvxy[j][i - 1] + mpw->pvxy[j - 1][i - 1]));
    float v = mpw->pvxy[j][i] + mpw->pvyx[j][i];

    /* The stress updates are stored in internal variables due to extensions necessary for FWI */
    /* calculate stress component update */
    u1 = ((mpm->pc55ipjp[j][i] * v) + (mpm->pc15ipjp[j][i] * vxxipjp) + (mpm->pc35ipjp[j][i] * vyyipjp));
    u2 = ((mpm->pc11[j][i] * mpw->pvxx[j][i]) + (mpm->pc13[j][i] * mpw->pvyy[j][i]) + (mpm->pc15[j][i] * vij));
    u3 = ((mpm->pc13[j][i] * mpw->pvxx[j][i]) + (mpm->pc33[j][i] * mpw->pvyy[j][i]) + (mpm->pc35[j][i] * vij));

    /* updating components of the stress tensor */
    mpw->psxy[j][i] += u1;
    mpw->psxx[j][i] += u2;
    mpw->psyy[j][i] += u3;
}
