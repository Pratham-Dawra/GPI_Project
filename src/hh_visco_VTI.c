
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
 *
 *   ------------------------------------------------------------- */

#include "fd.h"
#include "logging.h"

void model_visco_vti(float **rho, float **pc11, float **pc33, float **pc13, float **pc55,
                     float **ptau11, float **ptau33, float **ptau13, float **ptau55, float *eta, GlobVar *gv)
{
    float c11, c33, c55, c13, Rho, tau11, tau33, tau55, tau13;
    float sumc11, sumc13, sumc33, sumc55, ws, y;
    int ii, jj;
    char modfile[STRING_SIZE + 16];

    /*-----------------material property definition -------------------------*/

    /* anisotropic case, values from
     * Jones, Wang, 1981, Geophysics, 46, 3, 288-297 */

    const float C11 = 34.3e9, C33 = 22.7e9, C55 = 5.4e9, C13 = 10.7e9, RHO = 2000.0;
    const float TAU11 = gv->TAU, TAU33 = gv->TAU, TAU55 = gv->TAU, H = -1.0;

    /*-----------------------------------------------------------------------*/

    /* vector for maxwellbodies */
    float *pts = vector(1, gv->L);
    for (int l = 1; l <= gv->L; l++) {
        pts[l] = 1.0 / (2.0 * PI * gv->FL[l]);
        eta[l] = gv->DT / pts[l];
    }

    float fc = 1.0 / gv->TS;
    log_infoc(0, "VTI: center source frequency of %5.2fHz applied for calculation of relaxed moduli.\n", fc);

    ws = 2.0 * PI * fc;

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        for (int j = 1; j <= gv->NYG; j++) {

            y = (float)j *gv->DH;

            if (y < H) {
                c11 = 0.0;
                c33 = 0.0;
                c13 = 0.0;
                c55 = 1.0;
                Rho = 1.0;
                tau11 = 0.0;
                tau33 = 0.0;
                tau55 = 0.0;
            } else {
                c11 = C11;
                c33 = C33;
                c13 = C13;
                c55 = C55;
                Rho = RHO;
                tau11 = TAU11;
                tau33 = TAU33;
                tau55 = TAU55;
            }

            sumc11 = 0.0;
            sumc33 = 0.0;
            sumc55 = 0.0;
            for (int l = 1; l <= gv->L; l++) {
                sumc11 = sumc11 + ((ws * ws * pts[l] * pts[l] * tau11) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumc33 = sumc33 + ((ws * ws * pts[l] * pts[l] * tau33) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumc55 = sumc55 + ((ws * ws * pts[l] * pts[l] * tau55) / (1.0 + ws * ws * pts[l] * pts[l]));
            }

            /* relaxed moduli */
            c11 = c11 / (1.0 + sumc11);
            c33 = c33 / (1.0 + sumc33);
            c55 = c55 / (1.0 + sumc55);

            /* isotropic attenuation */
            tau13 = (c11 * gv->L * tau11 - 2.0 * c55 * gv->L * tau55) / (c11 - 2.0 * c55);

            sumc13 = 0.0;
            for (int l = 1; l <= gv->L; l++) {
                sumc13 = sumc13 + ((ws * ws * pts[l] * pts[l] * tau13) / (1.0 + ws * ws * pts[l] * pts[l]));
            }
            c13 = c13 / (1.0 + sumc13);

            /* only the PE which belongs to the current global gridpoint
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;

                pc11[jj][ii] = c11;
                rho[jj][ii] = Rho;
                pc33[jj][ii] = c33;
                pc13[jj][ii] = c13;
                pc55[jj][ii] = c55;

                ptau11[jj][ii] = tau11;
                ptau33[jj][ii] = tau33;
                ptau13[jj][ii] = tau13;
                ptau55[jj][ii] = tau55;
            }
        }
    }

    /* each PE writes his model to disk */

    /* all models are written to file */
    if (gv->WRITE_MODELFILES == 1) {
        sprintf(modfile, "%s.SOFI2D.c11", gv->MFILE);
        writemod(modfile, pc11, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c33", gv->MFILE);
        writemod(modfile, pc33, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c13", gv->MFILE);
        writemod(modfile, pc13, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c55", gv->MFILE);
        writemod(modfile, pc55, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau11", gv->MFILE);
        writemod(modfile, ptau11, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau33", gv->MFILE);
        writemod(modfile, ptau33, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau13", gv->MFILE);
        writemod(modfile, ptau13, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau55", gv->MFILE);
        writemod(modfile, ptau55, 3, gv);
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
}
