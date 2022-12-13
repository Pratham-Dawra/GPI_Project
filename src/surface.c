
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

/*------------------------------------------------------------------------
 *   stress free surface condition
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface(int ndepth, float **vx, float **vy, float **sxx, float **syy,
             float **sxy, float ***p, float ***q,
             float **pi, float **u, float **taup,
             float **taus, float *etajm, float *peta, float *hc, float *K_x, float *a_x, float *b_x, float **psi_vxx,
             GlobVar *gv)
{

    int h1;
    float bjm, djm, e, fjm, g;
    float vxx, vyy, sump = 0.0f;

    int fdoh = gv->FDORDER / 2;
    float dhi = 1.0f / gv->DH;
    float dthalbe = gv->DT / 2.0f;

    int j = ndepth;                 /* The free surface is located exactly in y=dh !! */
    for (int i = 1; i <= gv->NX; i++) {

        /* Compute values for shearmodulus u[j][i], P-wave modulus pi[j][i],
         * tau for S-waves and P-waves taus[j][i], taup[j][i] at half indizes: */

/*		for (l=1;l<=gv->L;l++){
			etajm[l]=0.5*(eta(i,j-1,l)+eta[j][i][l]);
		}
*/
        for (int l = 1; l <= gv->L; l++) {
            etajm[l] = peta[l];
        }

        /*Mirroring the components of the stress tensor to make
         * a stress free surface (method of imaging) */
        syy[j][i] = 0.0;

        /* since syy is zero on the free surface also the
         * corresponding memory-variables must set to zero */
        for (int l = 1; l <= gv->L; l++)
            q[j][i][l] = 0.0;

        /* now updating the stress component sxx and the memory-
         * variables p[j][i][l] at the free surface */

        /* first calculate spatial derivatives of components
         * of particle velocities */
        vxx = 0.0;
        vyy = 0.0;
        for (int m = 1; m <= fdoh; m++) {
            /*Mirroring the components of the stress tensor to make
             * a stress free surface (method of imaging) */
            syy[j - m][i] = -syy[j + m][i];
            sxy[j - m][i] = -sxy[j + m - 1][i];

            vxx += hc[m] * (vx[j][i + m - 1] - vx[j][i - m]);
            vyy += hc[m] * (vy[j + m - 1][i] - vy[j - m][i]);
        }
        vxx *= dhi;
        vyy *= dhi;

        if (gv->ABS_TYPE == 1) {
            /* apply PML boundary */
            /* left boundary */
            if ((!gv->BOUNDARY) && (gv->POS[1] == 0) && (i <= gv->FW)) {

                psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                vxx = vxx / K_x[i] + psi_vxx[j][i];
            }

            /* right boundary */
            if ((!gv->BOUNDARY) && (gv->POS[1] == gv->NPROCX - 1) && (i >= gv->NX - gv->FW + 1)) {

                h1 = (i - gv->NX + 2 * gv->FW);

                psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                vxx = vxx / K_x[h1] + psi_vxx[j][h1];
            }
        }

        /* sums used in updating sxx */
        sump = 0.0;
        for (int l = 1; l <= gv->L; l++)
            sump += p[j][i][l];

        fjm = u[j][i] * 2.0 * (1.0 + gv->L * taus[j][i]);
        g = pi[j][i] * (1.0 + gv->L * taup[j][i]);

        /* partially updating sxx */
        sxx[j][i] += -(gv->DT * (g - fjm) * (g - fjm) * vxx / g) - (gv->DT * (g - fjm) * vyy) - (dthalbe * sump);

        /* updating the memory-variable p[j][i][l] at the free surface */
        sump = 0.0;
        for (int l = 1; l <= gv->L; l++) {
            bjm = etajm[l] / (1.0 + (etajm[l] * 0.5));
            djm = 2.0 * u[j][i] * taus[j][i];
            e = pi[j][i] * taup[j][i];
            p[j][i][l] += bjm * (((djm - e) * ((fjm / g) - 1.0) * vxx) - ((djm - e) * vyy));
            sump += p[j][i][l];
        }
        /*completely updating the stress sxx */
        sxx[j][i] += (dthalbe * sump);
    }
}
