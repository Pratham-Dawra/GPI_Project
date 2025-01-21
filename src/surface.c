
/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
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
--------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * stress free surface condition
 *----------------------------------------------------------------------*/

#include "fd.h"

void surface(int ndepth, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv)
{
    int h1;
    float bjm, djm, e, fjm, g;
    float vxx, vyy, sump = 0.0f;

    int fdoh = gv->FDORDER / 2;
    float dhi = 1.0 / gv->DH;
    float dthalbe = gv->DT / 2.0;

    int j = ndepth;             /* The free surface is located exactly in y=0 */
    for (int i = 1; i <= gv->NX; i++) {

        /* Compute values for shearmodulus u[j][i], P-wave modulus pi[j][i],
         * tau for S-waves and P-waves taus[j][i], taup[j][i] at half indizes: */

        /* Mirroring the components of the stress tensor to make
         * a stress free surface (method of imaging) */
        mpw->psyy[j][i] = 0.0;

        /* since syy is zero on the free surface also the
         * corresponding memory-variables must set to zero */
        for (int l = 1; l <= gv->L; l++)
	    mpw->pq[j][i][l] = 0.0;

        /* now updating the stress component sxx and the memory-
         * variables p[j][i][l] at the free surface */

        /* first calculate spatial derivatives of components
         * of particle velocities */
        vxx = 0.0;
        vyy = 0.0;
        for (int m = 1; m <= fdoh; m++) {
            /* Mirroring the components of the stress tensor to make
             * a stress free surface (method of imaging) */
            mpw->psyy[j - m][i] = -mpw->psyy[j + m][i];
            mpw->psxy[j - m][i] = -mpw->psxy[j + m - 1][i];

            vxx += hc[m] * (mpw->pvx[j][i + m - 1] - mpw->pvx[j][i - m]);
            vyy += hc[m] * (mpw->pvy[j + m - 1][i] - mpw->pvy[j - m][i]);
        }
        vxx *= dhi;
        vyy *= dhi;

        if (gv->ABS_TYPE == 1) {
            /* apply PML boundary */
            /* left boundary */
            if ((!gv->BOUNDARY) && (gv->POS[1] == 0) && (i <= gv->FW)) {

                mpw->psi_vxx[j][i] = mpm->b_x[i] * mpw->psi_vxx[j][i] + mpm->a_x[i] * vxx;
                vxx = vxx / mpm->K_x[i] + mpw->psi_vxx[j][i];
            }

            /* right boundary */
            if ((!gv->BOUNDARY) && (gv->POS[1] == gv->NPROCX - 1) && (i >= gv->NX - gv->FW + 1)) {

                h1 = (i - gv->NX + 2 * gv->FW);

                mpw->psi_vxx[j][h1] = mpm->b_x[h1] * mpw->psi_vxx[j][h1] + mpm->a_x[h1] * vxx;
                vxx = vxx / mpm->K_x[h1] + mpw->psi_vxx[j][h1];
            }
        }

        /* sums used in updating sxx */
        sump = 0.0;
        for (int l = 1; l <= gv->L; l++)
            sump += mpw->pp[j][i][l];

        fjm = mpm->pu[j][i] * 2.0 * (1.0 + gv->L * mpm->ptaus[j][i]);
        g = mpm->ppi[j][i] * (1.0 + gv->L * mpm->ptaup[j][i]);

        /* partially updating sxx */
        mpw->psxx[j][i] += -(gv->DT * (g - fjm) * (g - fjm) * vxx / g) - (gv->DT * (g - fjm) * vyy) - (dthalbe * sump);

        /* updating the memory-variable p[j][i][l] at the free surface */
        sump = 0.0;
        for (int l = 1; l <= gv->L; l++) {
            bjm = mpm->peta[l] / (1.0 + (mpm->peta[l] * 0.5));
            djm = 2.0 * mpm->pu[j][i] * mpm->ptaus[j][i];
            e = mpm->ppi[j][i] * mpm->ptaup[j][i];
            mpw->pp[j][i][l] += bjm * (((djm - e) * ((fjm / g) - 1.0) * vxx) - ((djm - e) * vyy));
            sump += mpw->pp[j][i][l];
        }
        /* completely updating the stress sxx */
        mpw->psxx[j][i] += (dthalbe * sump);
    }
}
