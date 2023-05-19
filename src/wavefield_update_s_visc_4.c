
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

/* $Id: wavefield_update_s_visc.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the viscoelastic case*/

#include "fd.h"

void wavefield_update_s_visc_4(int i, int j, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;
    float sumr_old = 0.0f, sump_old = 0.0f, sumq_old = 0.0f;
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    float sumxx = 0.0, sumyy = 0.0, sumxy = 0.0, sumyx = 0.0;
    float u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;
    float ctemp;

    float dhi = 1.0 / gv->DH;
    float dthalbe = gv->DT / 2.0;

    // Save derviations
    mpw->vxx_1[j][i] = mpw->pvxx[j][i] * gv->DH;
    mpw->vyy_1[j][i] = mpw->pvyy[j][i] * gv->DH;
    mpw->vxy_1[j][i] = mpw->pvxy[j][i] * gv->DH;
    mpw->vyx_1[j][i] = mpw->pvyx[j][i] * gv->DH;

    // Calculate Adams-Bashforth stuff
    sumxx = c1 * mpw->vxx_1[j][i] + c2 * mpw->vxx_2[j][i] + c3 * mpw->vxx_3[j][i] + c4 * mpw->vxx_4[j][i];
    sumyy = c1 * mpw->vyy_1[j][i] + c2 * mpw->vyy_2[j][i] + c3 * mpw->vyy_3[j][i] + c4 * mpw->vyy_4[j][i];
    sumxy = c1 * mpw->vxy_1[j][i] + c2 * mpw->vxy_2[j][i] + c3 * mpw->vxy_3[j][i] + c4 * mpw->vxy_4[j][i];
    sumyx = c1 * mpw->vyx_1[j][i] + c2 * mpw->vyx_2[j][i] + c3 * mpw->vyx_3[j][i] + c4 * mpw->vyx_4[j][i];

    for (int l = 1; l <= gv->L; l++) {
        /* computing sums of the old memory variables */
        sumr_old += c1 * mpw->pr[j][i][l] + c2 * mpw->pr_2[j][i][l] + c3 * mpw->pr_3[j][i][l] + c4 * mpw->pr_4[j][i][l];
        sump_old += c1 * mpw->pp[j][i][l] + c2 * mpw->pp_2[j][i][l] + c3 * mpw->pp_3[j][i][l] + c4 * mpw->pp_4[j][i][l];
        sumq_old += c1 * mpw->pq[j][i][l] + c2 * mpw->pq_2[j][i][l] + c3 * mpw->pq_3[j][i][l] + c4 * mpw->pq_4[j][i][l];
        /* now updating the memory-variables and sum them up */
        ctemp = 2 * (1 - mpm->cip[l]) / c1;
        mpw->pr_4[j][i][l] =
            mpm->bip[l] * (mpw->pr[j][i][l] * mpm->cip[l] - ctemp * c2 * (mpw->pr[j][i][l] + mpw->pr_2[j][i][l]) -
                           ctemp * c3 * (mpw->pr_2[j][i][l] + mpw->pr_3[j][i][l]) - ctemp * c4 * (mpw->pr_3[j][i][l] +
                                                                                                  mpw->pr_4[j][i][l]) -
                           (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
        mpw->pp_4[j][i][l] =
            mpm->bjm[l] * (mpw->pp[j][i][l] * mpm->cjm[l] - ctemp * c2 * (mpw->pp[j][i][l] + mpw->pp_2[j][i][l]) -
                           ctemp * c3 * (mpw->pp_2[j][i][l] + mpw->pp_3[j][i][l]) - ctemp * c4 * (mpw->pp_3[j][i][l] +
                                                                                                  mpw->pp_4[j][i][l]) -
                           (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) + (2.0 * mpm->d[j][i][l] * sumyy * dhi));
        mpw->pq_4[j][i][l] =
            mpm->bjm[l] * (mpw->pq[j][i][l] * mpm->cjm[l] - ctemp * c2 * (mpw->pq[j][i][l] + mpw->pq_2[j][i][l]) -
                           ctemp * c3 * (mpw->pq_2[j][i][l] + mpw->pq_3[j][i][l]) - ctemp * c4 * (mpw->pq_3[j][i][l] +
                                                                                                  mpw->pq_4[j][i][l]) -
                           (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) + (2.0 * mpm->d[j][i][l] * sumxx * dhi));
        sumr += c1 * mpw->pr_4[j][i][l] + c2 * mpw->pr[j][i][l] + c3 * mpw->pr_2[j][i][l] + c4 * mpw->pr_3[j][i][l];
        sump += c1 * mpw->pp_4[j][i][l] + c2 * mpw->pp[j][i][l] + c3 * mpw->pp_2[j][i][l] + c4 * mpw->pp_3[j][i][l];
        sumq += c1 * mpw->pq_4[j][i][l] + c2 * mpw->pq[j][i][l] + c3 * mpw->pq_2[j][i][l] + c4 * mpw->pq_3[j][i][l];
    }

    /* The stress updates are stored in internal variables due to extensions necessary for FWI */
    /* calculate stress component update */
    u1 = (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * (sumr_old + sumr));
    u2 = (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * (sump_old + sump));
    u3 = (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * (sumq_old + sumq));

    /* updating components of the stress tensor */
    mpw->psxy[j][i] += u1;
    mpw->psxx[j][i] += u2;
    mpw->psyy[j][i] += u3;
}
