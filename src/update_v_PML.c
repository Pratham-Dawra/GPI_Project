
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_v_PML(int nx2, int ny2, int nt, int sw, MemModel *mpm, MemWavefield *mpw, MemInv * minv,
                  GlobVar *gv, GlobVarInv *vinv)
{
    int h1;
    float sxx_x, syy_y, sxy_y, sxy_x;
    double time1 = 0.0, time2 = 0.0;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating particle velocities...\n");
    }

    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */

  /*-------------------------update-CPML-Boundarys------------------------------------*/

    /* left boundary */
    for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            cpml_update_v_x(i, j, &sxx_x, &sxy_x, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* right boundary */
    for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_v_x(h1, j, &sxx_x, &sxy_x, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* top boundary */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            cpml_update_v_y(i, j, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* bottom boundary */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_v_y(i, h1, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* corners */

    /*left-top */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            cpml_update_v_x(i, j, &sxx_x, &sxy_x, mpm, mpw);
            cpml_update_v_y(i, j, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /*left-bottom */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            cpml_update_v_x(i, j, &sxx_x, &sxy_x, mpm, mpw);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_v_y(i, h1, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* right-top */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_v_x(h1, j, &sxx_x, &sxy_x, mpm, mpw);
            cpml_update_v_y(i, j, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    /* right-bottom */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_v_x(h1, j, &sxx_x, &sxy_x, mpm, mpw);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_v_y(i, h1, &sxy_y, &syy_y, mpm, mpw);
            wavefield_update_v(i, j, sw, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, minv, gv, vinv);
        }
    }

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating particle velocities (real time: %4.3fs).\n", time2 - time1);
    }
}
