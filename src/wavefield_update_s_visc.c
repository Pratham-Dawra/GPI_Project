
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

void wavefield_update_s_visc(int i, int j, MemModel * mpm, MemWavefield * mpw, MemInv *minv, GlobVar * gv)
{
    float dthalbe = gv->DT / 2.0;
    float sumr_old = 0.0f, sump_old = 0.0f, sumq_old = 0.0f;
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;

    for (int l = 1; l <= gv->L; l++) {
        /* computing sums of the old memory variables */
        sumr_old += mpw->pr[j][i][l];
        sump_old += mpw->pp[j][i][l];
        sumq_old += mpw->pq[j][i][l];
        /* now updating the memory-variables and sum them up */
        mpw->pr[j][i][l] = mpm->bip[l] * (mpw->pr[j][i][l] * mpm->cip[l] - (mpm->dip[j][i][l] * (mpw->pvxy[j][i] + mpw->pvyx[j][i])));
        mpw->pp[j][i][l] =
            mpm->bjm[l] * (mpw->pp[j][i][l] * mpm->cjm[l] - (mpm->e[j][i][l] * (mpw->pvxx[j][i] + mpw->pvyy[j][i])) +
                           (2.0 * mpm->d[j][i][l] * mpw->pvyy[j][i]));
        mpw->pq[j][i][l] =
            mpm->bjm[l] * (mpw->pq[j][i][l] * mpm->cjm[l] - (mpm->e[j][i][l] * (mpw->pvxx[j][i] + mpw->pvyy[j][i])) +
                           (2.0 * mpm->d[j][i][l] * mpw->pvxx[j][i]));
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    /* calculate stress component update */
    u1 = (mpm->fipjp[j][i] * (mpw->pvxy[j][i] + mpw->pvyx[j][i])) + (dthalbe * (sumr_old + sumr));
    u2 = (mpm->g[j][i] * (mpw->pvxx[j][i] + mpw->pvyy[j][i])) - (2.0 * mpm->f[j][i] * mpw->pvyy[j][i]) + (dthalbe * (sump_old + sump));
    u3 = (mpm->g[j][i] * (mpw->pvxx[j][i] + mpw->pvyy[j][i])) - (2.0 * mpm->f[j][i] * mpw->pvxx[j][i]) + (dthalbe * (sumq_old + sumq));

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
