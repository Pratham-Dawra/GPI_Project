
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

    float vxxipjp, vyyipjp, vyxij, vxyij, vij, v;

    /* computing sums of the old memory variables */
    float dthalbe = gv->DT / 2.0;
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    for (int l = 1; l <= gv->L; l++) {
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    vxxipjp = 0.25 * (mpw->pvxx[j][i] + mpw->pvxx[j + 1][i] + mpw->pvxx[j][i + 1] + mpw->pvxx[j + 1][i + 1]);
    vyyipjp = 0.25 * (mpw->pvyy[j][i] + mpw->pvyy[j + 1][i] + mpw->pvyy[j][i + 1] + mpw->pvyy[j + 1][i + 1]);

    vyxij = 0.25 * (mpw->pvyx[j][i] + mpw->pvyx[j - 1][i] + mpw->pvyx[j][i - 1] + mpw->pvyx[j - 1][i - 1]);
    vxyij = 0.25 * (mpw->pvxy[j][i] + mpw->pvxy[j - 1][i] + mpw->pvxy[j][i - 1] + mpw->pvxy[j - 1][i - 1]);
    vij = vyxij + vxyij;

    v = mpw->pvxy[j][i] + mpw->pvyx[j][i];

    /* updating components of the stress tensor, partially */
    mpw->psxy[j][i] +=
        (mpm->pc55ipjpu[j][i] * v) + (mpm->pc15ipjpu[j][i] * vxxipjp) + (mpm->pc35ipjpu[j][i] * vyyipjp) +
        (dthalbe * sumr);
    mpw->psxx[j][i] +=
        (mpm->pc11u[j][i] * mpw->pvxx[j][i]) + (mpm->pc13u[j][i] * mpw->pvyy[j][i]) + +(mpm->pc15u[j][i] * vij) +
        (dthalbe * sump);
    mpw->psyy[j][i] +=
        (mpm->pc13u[j][i] * mpw->pvxx[j][i]) + (mpm->pc33u[j][i] * mpw->pvyy[j][i]) + (mpm->pc35u[j][i] * vij) +
        (dthalbe * sumq);

    /* now updating the memory-variables and sum them up */
    sumr = sump = sumq = 0.0f;
    for (int l = 1; l <= gv->L; l++) {
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

    /* and now the components of the stress tensor are completely updated */
    mpw->psxy[j][i] += (dthalbe * sumr);
    mpw->psxx[j][i] += (dthalbe * sump);
    mpw->psyy[j][i] += (dthalbe * sumq);

}
