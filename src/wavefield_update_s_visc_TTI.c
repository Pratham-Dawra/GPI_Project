
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
 * Update Function of the stress-Wavefields in the viscoelastic case
 */

#include "fd.h"

void wavefield_update_s_visc_TTI(int i, int j, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{

    float dthalbe = gv->DT / 2.0;
    float sumr_old = 0.0f, sump_old = 0.0f, sumq_old = 0.0f;
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    float vxxipjp, vyyipjp, vyxij, vxyij, vij, v;
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;

    vxxipjp = 0.25 * (mpw->pvxx[j][i] + mpw->pvxx[j + 1][i] + mpw->pvxx[j][i + 1] + mpw->pvxx[j + 1][i + 1]);
    vyyipjp = 0.25 * (mpw->pvyy[j][i] + mpw->pvyy[j + 1][i] + mpw->pvyy[j][i + 1] + mpw->pvyy[j + 1][i + 1]);

    vyxij = 0.25 * (mpw->pvyx[j][i] + mpw->pvyx[j - 1][i] + mpw->pvyx[j][i - 1] + mpw->pvyx[j - 1][i - 1]);
    vxyij = 0.25 * (mpw->pvxy[j][i] + mpw->pvxy[j - 1][i] + mpw->pvxy[j][i - 1] + mpw->pvxy[j - 1][i - 1]);
    vij = vyxij + vxyij;

    v = mpw->pvxy[j][i] + mpw->pvyx[j][i];

    for (int l = 1; l <= gv->L; l++) {
        /* computing sums of the old memory variables */
        sumr_old += mpw->pr[j][i][l];
        sump_old += mpw->pp[j][i][l];
        sumq_old += mpw->pq[j][i][l];
        /* now updating the memory-variables and sum them up */
        mpw->pr[j][i][l] =
            mpm->bip[l] * (mpw->pr[j][i][l] * mpm->cip[l] - (mpm->pc55ipjpd[j][i][l] * v) -
                           (mpm->pc15ipjpd[j][i][l] * vxxipjp) - (mpm->pc35ipjpd[j][i][l] * vyyipjp));
        mpw->pp[j][i][l] =
            mpm->bip[l] * (mpw->pp[j][i][l] * mpm->cip[l] - (mpm->pc11d[j][i][l] * mpw->pvxx[j][i]) -
                           (mpm->pc13d[j][i][l] * mpw->pvyy[j][i]) - (mpm->pc15d[j][i][l] * vij));
        mpw->pq[j][i][l] =
            mpm->bip[l] * (mpw->pq[j][i][l] * mpm->cip[l] - (mpm->pc13d[j][i][l] * mpw->pvxx[j][i]) -
                           (mpm->pc33d[j][i][l] * mpw->pvyy[j][i]) - (mpm->pc35d[j][i][l] * vij));
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    /* The stress updates are stored in internal variables due to extensions necessary for FWI */
    /* calculate stress component update */
    u1 = (mpm->pc55ipjpu[j][i] * v) + (mpm->pc15ipjpu[j][i] * vxxipjp) + (mpm->pc35ipjpu[j][i] * vyyipjp) +
    (dthalbe * (sumr_old + sumr));
    u2 = (mpm->pc11u[j][i] * mpw->pvxx[j][i]) + (mpm->pc13u[j][i] * mpw->pvyy[j][i]) + +(mpm->pc15u[j][i] * vij) +
    (dthalbe * (sump_old + sump));
    u3 = (mpm->pc13u[j][i] * mpw->pvxx[j][i]) + (mpm->pc33u[j][i] * mpw->pvyy[j][i]) + (mpm->pc35u[j][i] * vij) +
    (dthalbe * (sumq_old + sumq));

    /* updating components of the stress tensor */
    mpw->psxy[j][i] += u1;
    mpw->psxx[j][i] += u2;
    mpw->psyy[j][i] += u3;
}
