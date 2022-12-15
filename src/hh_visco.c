
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

/* -------------------------------------------------------------
 *   Model homogeneous half space
 *   ------------------------------------------------------------- */

#include "fd.h"
#include "logging.h"

void model_visco(float **rho, float **pi, float **u, float **taus, float **taup, float *eta, GlobVar *gv)
{
    float Rhov, muv, piv, Vp, Vs, y;
    float ts, tp, sumu, sumpi;
    int ii, jj;
    char modfile[STRING_SIZE + 16];
    float **pwavemod = NULL, **swavemod = NULL;

    /*-----------------material property definition -------------------------*/

    /* parameters for layer 1 */
    const float vp1 = 0.0, vs1 = 1.0, rho1 = 1.0, tp1 = 0.01, ts1 = 0.01, h = 100.0;

    /* parameters for layer 2 */
    const float vp2 = 3500.0, vs2 = 2000.0, rho2 = 2000.0, tp2 = 0.01, ts2 = 0.01;;
    const float h1 = 95.0, h2 = 110.0, tlp = 2.0, tls = 2.0;
    /*const float h1=95.0, h2=110.0, tlp=2.0, tls=2.0; */

    /*-----------------------------------------------------------------------*/

    if (gv->WRITE_MODELFILES == 1) {
        pwavemod = matrix(0, gv->NY + 1, 0, gv->NX + 1);
        swavemod = matrix(0, gv->NY + 1, 0, gv->NX + 1);
    }

    /* vector for maxwellbodies */
    float *pts = vector(1, gv->L);
    for (int l = 1; l <= gv->L; l++) {
        pts[l] = 1.0 / (2.0 * PI * gv->FL[l]);
        eta[l] = gv->DT / pts[l];
    }

    float fc = 1.0 / gv->TS;
    log_infoc(0, "Center source frequency of %5.2fHz applied for calculation of relaxed moduli.\n", fc);

    float ws = 2.0 * PI * fc;

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        for (int j = 1; j <= gv->NYG; j++) {

            /* calculate coordinate in m */
            y = (float)j *gv->DH;

            /* two layer case */
            if (y <= h) {
                Vp = vp1;
                Vs = vs1;
                Rhov = rho1;
                tp = tp1;
                ts = ts1;
            }

            else {
                Vp = vp2;
                Vs = vs2;
                Rhov = rho2;
                tp = tp2;
                ts = ts2;
            }

            if ((y >= h1) && (y <= h2)) {
                tp = tlp;
                ts = tls;
            }

            sumu = 0.0;
            sumpi = 0.0;
            for (int l = 1; l <= gv->L; l++) {
                sumu = sumu + ((ws * ws * pts[l] * pts[l] * ts) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumpi = sumpi + ((ws * ws * pts[l] * pts[l] * tp) / (1.0 + ws * ws * pts[l] * pts[l]));
            }

            muv = Vs * Vs * Rhov / (1.0 + sumu);
            piv = Vp * Vp * Rhov / (1.0 + sumpi);

            /* only the PE which belongs to the current global gridpoint
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;

                taus[jj][ii] = ts;
                taup[jj][ii] = tp;
                u[jj][ii] = muv;
                rho[jj][ii] = Rhov;
                pi[jj][ii] = piv;
                if (gv->WRITE_MODELFILES == 1) {
                    pwavemod[jj][ii] = Vp;
                    swavemod[jj][ii] = Vs;
                }
            }
        }
    }

    /* each PE writes his model to disk */

    /* all models are written to file */
    if (gv->WRITE_MODELFILES == 1) {
        sprintf(modfile, "%s.SOFI2D.u", gv->MFILE);
        writemod(modfile, u, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.pi", gv->MFILE);
        writemod(modfile, pi, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.ts", gv->MFILE);
        writemod(modfile, taus, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tp", gv->MFILE);
        writemod(modfile, taup, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.vp", gv->MFILE);
        writemod(modfile, pwavemod, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.vs", gv->MFILE);
        writemod(modfile, swavemod, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, rho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }

    /* only density is written to file */
    if (gv->WRITE_MODELFILES == 2) {
        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, rho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }

    free_vector(pts, 1, gv->L);
    if (gv->WRITE_MODELFILES == 1) {
        free_matrix(pwavemod, 0, gv->NY + 1, 0, gv->NX + 1);
        free_matrix(swavemod, 0, gv->NY + 1, 0, gv->NX + 1);
    }
}
