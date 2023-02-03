
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

void model_visco_tti(MemModel * mpm, GlobVar * gv)
{
    float c11, c33, c55, c13, Rho;
    float tau11, tau13, tau33, tau55, tau15, tau35;
    float sumc11, sumc13, sumc33, sumc55;
    float a1, a3, a4, a5, a6;
    float c11t, c33t, c55t, c13t, c15t, c35t;
    int ii, jj;
    char modfile[STRING_SIZE + 16];

    /*-----------------material property definition -------------------------*/

    /* anisotropic case, values from
     * Jones, Wang, 1981, Geophysics, 46, 3, 288-297 */

    const float C11 = 34.3e9, C33 = 22.7e9, C55 = 5.4e9, C13 = 10.7e9, RHO = 2000.0;
    const float TAU11 = gv->TAU, TAU33 = gv->TAU, TAU55 = 1.0 * gv->TAU;
    const float THETA = 30.0;   /* rotation angle in degrees */
    const int jh = -1;

    /*-----------------------------------------------------------------------*/

    float t = THETA * PI / 180.0;

    float l1 = cos(t);
    float l2 = sin(t);
    float l12 = l1 * l1;
    float l22 = l2 * l2;
    float l14 = l12 * l12;
    float l24 = l22 * l22;
    float l13 = l1 * l12;
    float l23 = l2 * l22;

    /* vector for maxwellbodies */
    float *pts = vector(1, gv->L);
    for (int l = 1; l <= gv->L; l++) {
        pts[l] = 1.0 / (2.0 * PI * gv->FL[l]);
        mpm->peta[l] = gv->DT / pts[l];
    }

    //float fc = 1.0 / gv->TS;
    float ws = 2.0 * PI * gv->F_REF;
    log_infoc(0, "TTI: Center frequency of %5.2fHz applied for calculation of relaxed moduli.\n", gv->F_REF);

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        for (int j = 1; j <= gv->NYG; j++) {

            if (j < jh) {
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
            tau13 = tau15 = tau35 = (c11 * gv->L * tau11 - 2.0 * c55 * gv->L * tau55) / (c11 - 2.0 * c55);

            sumc13 = 0.0;
            for (int l = 1; l <= gv->L; l++) {
                sumc13 = sumc13 + ((ws * ws * pts[l] * pts[l] * tau13) / (1.0 + ws * ws * pts[l] * pts[l]));
            }
            c13 = c13 / (1.0 + sumc13);

            /* Bond transformation (Oh et al, 2020, GJI, doi: 10.1093/gji/ggaa295 */

            a1 = 2.0 * c13 + 4.0 * c55;
            a3 = c11 + c33 - 4.0 * c55;
            a4 = c11 + c33 - 2.0 * c13;
            a5 = c13 - c11 + 2.0 * c55;
            a6 = c13 - c33 + 2.0 * c55;

            c11t = c11 * l14 + c33 * l24 + a1 * l12 * l22;
            c33t = c11 * l24 + c33 * l14 + a1 * l12 * l22;
            c13t = a3 * l12 * l22 + c13 * (l14 + l24);
            c55t = a4 * l12 * l22 + c55 * (l12 - l22) * (l12 - l22);
            c15t = a5 * l13 * l2 - a6 * l1 * l23;
            c35t = a5 * l23 * l1 - a6 * l2 * l13;

            /* only the PE which belongs to the current global gridpoint
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;

                mpm->pc11[jj][ii] = c11t;
                mpm->prho[jj][ii] = Rho;
                mpm->pc33[jj][ii] = c33t;
                mpm->pc13[jj][ii] = c13t;
                mpm->pc55[jj][ii] = c55t;
                mpm->pc15[jj][ii] = c15t;
                mpm->pc35[jj][ii] = c35t;

                mpm->ptau11[jj][ii] = tau11;
                mpm->ptau33[jj][ii] = tau33;
                mpm->ptau13[jj][ii] = tau13;
                mpm->ptau55[jj][ii] = tau55;
                mpm->ptau15[jj][ii] = tau15;
                mpm->ptau35[jj][ii] = tau35;

            }
        }
    }

    /* each PE writes his model to disk */

    /* all models are written to file */
    if (gv->WRITE_MODELFILES == 1) {
        sprintf(modfile, "%s.SOFI2D.c11", gv->MFILE);
        writemod(modfile, mpm->pc11, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c33", gv->MFILE);
        writemod(modfile, mpm->pc33, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c13", gv->MFILE);
        writemod(modfile, mpm->pc13, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c55", gv->MFILE);
        writemod(modfile, mpm->pc55, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c15", gv->MFILE);
        writemod(modfile, mpm->pc15, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.c35", gv->MFILE);
        writemod(modfile, mpm->pc35, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau11", gv->MFILE);
        writemod(modfile, mpm->ptau11, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau33", gv->MFILE);
        writemod(modfile, mpm->ptau33, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau13", gv->MFILE);
        writemod(modfile, mpm->ptau13, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau55", gv->MFILE);
        writemod(modfile, mpm->ptau55, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau15", gv->MFILE);
        writemod(modfile, mpm->ptau15, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.tau35", gv->MFILE);
        writemod(modfile, mpm->ptau35, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);

        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, mpm->prho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }

    /* only density is written to file */
    if (gv->WRITE_MODELFILES == 2) {
        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, mpm->prho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }
    free_vector(pts, 1, gv->L);
}
