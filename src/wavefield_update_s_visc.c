
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

/* (int i, int j, float vxx, float vyx, float vxy, float vyy, float **sxy,
                             float **sxx, float **syy, float ***r, float ***p,
                             float ***q, float **fipjp, float **f, float **g, float *bip,
                             float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, GlobVar * gv) */

void wavefield_update_s_visc(int i, int j, float vxx, float vyx, float vxy, float vyy, MemModel * mpm,
                             MemWavefield * mpw, MemInv *minv, GlobVar * gv)
{
    /* computing sums of the old memory variables */
    float dthalbe = gv->DT / 2.0;
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    for (int l = 1; l <= gv->L; l++) {
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    /* updating components of the stress tensor, partially */
    mpw->psxy[j][i] += (mpm->fipjp[j][i] * (vxy + vyx)) + (dthalbe * sumr);
    mpw->psxx[j][i] += (mpm->g[j][i] * (vxx + vyy)) - (2.0 * mpm->f[j][i] * vyy) + (dthalbe * sump);
    mpw->psyy[j][i] += (mpm->g[j][i] * (vxx + vyy)) - (2.0 * mpm->f[j][i] * vxx) + (dthalbe * sumq);

    /*if (gv->MODE == FWI) {
        minv->uxy[j][i] = mpw->psxy[j][i] / gv->DT;
        minv->ux[j][i] = mpw->psxx[j][i] / gv->DT;
        minv->uy[j][i] = mpw->psyy[j][i] / gv->DT;
    }*/

    /* now updating the memory-variables and sum them up */
    sumr = sump = sumq = 0.0f;
    for (int l = 1; l <= gv->L; l++) {
        mpw->pr[j][i][l] = mpm->bip[l] * (mpw->pr[j][i][l] * mpm->cip[l] - (mpm->dip[j][i][l] * (vxy + vyx)));
        mpw->pp[j][i][l] =
            mpm->bjm[l] * (mpw->pp[j][i][l] * mpm->cjm[l] - (mpm->e[j][i][l] * (vxx + vyy)) +
                           (2.0 * mpm->d[j][i][l] * vyy));
        mpw->pq[j][i][l] =
            mpm->bjm[l] * (mpw->pq[j][i][l] * mpm->cjm[l] - (mpm->e[j][i][l] * (vxx + vyy)) +
                           (2.0 * mpm->d[j][i][l] * vxx));
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    /* and now the components of the stress tensor are completely updated */
    mpw->psxy[j][i] += (dthalbe * sumr);
    mpw->psxx[j][i] += (dthalbe * sump);
    mpw->psyy[j][i] += (dthalbe * sumq);

    if (gv->MODE == FWI) {
        minv->uxy[j][i] = mpw->psxy[j][i] / gv->DT;
        minv->ux[j][i] = mpw->psxx[j][i] / gv->DT;
        minv->uy[j][i] = mpw->psyy[j][i] / gv->DT;
/*        minv->uxy[j][i] += (0.5 * sumr);
        minv->ux[j][i] += (0.5 * sump);
        minv->uy[j][i] += (0.5 * sumq); */
    }
}
