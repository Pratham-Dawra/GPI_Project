
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

/* void wavefield_update_s_visc_VTI(int i, int j, float vxx, float vyx, float vxy, float vyy,
                                 float **sxy, float **sxx, float **syy, float ***p, float ***r, float ***q,
                                 float **pc55ipjpu, float **pc13u, float **pc11u, float **pc33u,
                                 float ***pc55ipjpd, float ***pc13d, float ***pc11d, float ***pc33d,
                                 float *bip, float *cip, GlobVar *gv) */

void wavefield_update_s_visc_VTI(int i, int j, float vxx, float vyx, float vxy, float vyy, MemModel * mpm,
                                 MemWavefield * mpw, GlobVar * gv)
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
    mpw->psxy[j][i] += (mpm->pc55ipjpu[j][i] * (vxy + vyx)) + (dthalbe * sumr);
    mpw->psxx[j][i] += (mpm->pc11u[j][i] * vxx) + (mpm->pc13u[j][i] * vyy) + (dthalbe * sump);
    mpw->psyy[j][i] += (mpm->pc13u[j][i] * vxx) + (mpm->pc33u[j][i] * vyy) + (dthalbe * sumq);

    /* now updating the memory-variables and sum them up */
    sumr = sump = sumq = 0.0f;
    for (int l = 1; l <= gv->L; l++) {
        mpw->pr[j][i][l] = mpm->bip[l] * (mpw->pr[j][i][l] * mpm->cip[l] - (mpm->pc55ipjpd[j][i][l] * (vxy + vyx)));
        mpw->pp[j][i][l] =
            mpm->bip[l] * (mpw->pp[j][i][l] * mpm->cip[l] - (mpm->pc11d[j][i][l] * vxx) - (mpm->pc13d[j][i][l] * vyy));
        mpw->pq[j][i][l] =
            mpm->bip[l] * (mpw->pq[j][i][l] * mpm->cip[l] - (mpm->pc13d[j][i][l] * vxx) - (mpm->pc33d[j][i][l] * vyy));
        sumr += mpw->pr[j][i][l];
        sump += mpw->pp[j][i][l];
        sumq += mpw->pq[j][i][l];
    }

    /* and now the components of the stress tensor are
     * completely updated */
    mpw->psxy[j][i] += (dthalbe * sumr);
    mpw->psxx[j][i] += (dthalbe * sump);
    mpw->psyy[j][i] += (dthalbe * sumq);

}
