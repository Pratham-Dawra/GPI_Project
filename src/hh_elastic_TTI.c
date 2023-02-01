
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
 *   if variable "h" is decreased, a layer over half-space is gained
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_elastic_TTI(MemModel * mpm, GlobVar * gv)
{
    float c11, c33, c55, c13, Rho;
    int ii, jj;
    float a1, a3, a4, a5, a6;
    float c11t, c33t, c55t, c13t, c15t, c35t;

    char modfile[STRING_SIZE + 16];

    /*-----------------material property definition -------------------------*/

    /* anisotropic case, values from
     * Jones, Wang, 1981, Geophysics, 46, 3, 288-297 */

    const float C11 = 34.3e9, C33 = 22.7e9, C55 = 5.4e9, C13 = 10.7e9, RHO = 2000.0;
    const float THETA = 0.0;    /* rotation angle in degrees */
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

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        for (int j = 1; j <= gv->NYG; j++) {

            if (j < jh) {
                c11 = 0.0;
                c33 = 0.0;
                c13 = 0.0;
                c55 = 1.0;
                Rho = 1.0;
            } else {
                c11 = C11;
                c33 = C33;
                c13 = C13;
                c55 = C55;
                Rho = RHO;
            }

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

            }
        }
    }

    /* each PE writes his model to disk */

    /* only the density model is written to file */
    if (gv->WRITE_MODELFILES == 2) {
        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, mpm->prho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }

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

        sprintf(modfile, "%s.SOFI2D.rho", gv->MFILE);
        writemod(modfile, mpm->prho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID == 0)
            mergemod(modfile, 3, gv);
    }
}
