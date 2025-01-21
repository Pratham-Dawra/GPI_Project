
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

void surface_elastic(int ndepth, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv)
{

    int h1;
    float fjm, g;
    float vxx, vyy;

    int fdoh = gv->FDORDER / 2;
    float dhi = 1.0 / gv->DH;

    int j = ndepth;             /* The free surface is located exactly in y=0 */
    for (int i = 1; i <= gv->NX; i++) {

        /* Mirroring the components of the stress tensor to make
         * a stress free surface (method of imaging) */
        mpw->psyy[j][i] = 0.0;

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

        fjm = mpm->pu[j][i] * 2.0;
        g = mpm->ppi[j][i];

        /* Update sxx without vertical derivates (last update will be canceld)
         * sxx=sxx_new - sxx =  gv->DT*fjm*(2-fjm/g)vxx   -  ( g* ( vxx+vyy ) - fjm *vyy)
         *                   = -(gv->DT*(g-fmj)*(g-fmj)*vxx/g)-(gv->DT*(g-fjm)*vyy) */

        mpw->psxx[j][i] += -(gv->DT * (g - fjm) * (g - fjm) * vxx / g) - (gv->DT * (g - fjm) * vyy);

    }
}
