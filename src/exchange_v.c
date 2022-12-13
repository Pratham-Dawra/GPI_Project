
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
 *   write values of dynamic field variables at the edges of the
 *   local grid into buffer arrays and  exchanged between
 *   processes.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_v(int nd, float **vx, float **vy,
                float **bufferlef_to_rig, float **bufferrig_to_lef,
                float **buffertop_to_bot, float **bufferbot_to_top, GlobVar *gv)
{

    MPI_Status status;
    int n;
    int fdo = nd;
    int fdo2 = 2 * fdo;

    /* top - bottom */

    if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
        for (int i = 1; i <= gv->NX; i++) {
            n = 1;
            /* storage of top of local volume into buffer */
            for (int l = 1; l <= fdo - 1; l++) {
                buffertop_to_bot[i][n++] = vy[l][i];
            }
            for (int l = 1; l <= fdo; l++) {
                buffertop_to_bot[i][n++] = vx[l][i];
            }
        }

    if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
        for (int i = 1; i <= gv->NX; i++) {
            /* storage of bottom of local volume into buffer */
            n = 1;
            for (int l = 1; l <= fdo; l++) {
                bufferbot_to_top[i][n++] = vy[gv->NY - l + 1][i];
            }
            for (int l = 1; l <= fdo - 1; l++) {
                bufferbot_to_top[i][n++] = vx[gv->NY - l + 1][i];
            }
        }

    /* send and reveive values for points at inner boundaries */
    MPI_Sendrecv_replace(&buffertop_to_bot[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[3], gv->TAG5, gv->INDEX[4],
                         gv->TAG5, MPI_COMM_WORLD, &status);
    MPI_Sendrecv_replace(&bufferbot_to_top[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[4], gv->TAG6, gv->INDEX[3],
                         gv->TAG6, MPI_COMM_WORLD, &status);

    if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
        for (int i = 1; i <= gv->NX; i++) {
            n = 1;
            for (int l = 1; l <= fdo - 1; l++) {
                vy[gv->NY + l][i] = buffertop_to_bot[i][n++];
            }
            for (int l = 1; l <= fdo; l++) {
                vx[gv->NY + l][i] = buffertop_to_bot[i][n++];
            }
        }

    if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
        for (int i = 1; i <= gv->NX; i++) {
            n = 1;
            for (int l = 1; l <= fdo; l++) {
                vy[1 - l][i] = bufferbot_to_top[i][n++];
            }
            for (int l = 1; l <= fdo - 1; l++) {
                vx[1 - l][i] = bufferbot_to_top[i][n++];
            }
        }

    /* left - right */

    /* exchange if periodic boundary condition is applied */
    if ((gv->BOUNDARY) || (gv->POS[1] != 0))
        for (int j = 1; j <= gv->NY; j++) {
            /* storage of left edge of local volume into buffer */
            n = 1;
            for (int l = 1; l <= fdo; l++) {
                bufferlef_to_rig[j][n++] = vy[j][l];
            }
            for (int l = 1; l <= fdo - 1; l++) {
                bufferlef_to_rig[j][n++] = vx[j][l];
            }
        }

    /* no exchange if periodic boundary condition is applied */
    if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
        for (int j = 1; j <= gv->NY; j++) {
            /* storage of right edge of local volume into buffer */
            n = 1;
            for (int l = 1; l <= fdo - 1; l++) {
                bufferrig_to_lef[j][n++] = vy[j][gv->NX - l + 1];
            }
            for (int l = 1; l <= fdo; l++) {
                bufferrig_to_lef[j][n++] = vx[j][gv->NX - l + 1];
            }
        }

    /* send and reveive values for points at inner boundaries */
    MPI_Sendrecv_replace(&bufferlef_to_rig[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[1], gv->TAG1, gv->INDEX[2],
                         gv->TAG1, MPI_COMM_WORLD, &status);
    MPI_Sendrecv_replace(&bufferrig_to_lef[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[2], gv->TAG2, gv->INDEX[1],
                         gv->TAG2, MPI_COMM_WORLD, &status);

    /* no exchange if periodic boundary condition is applied */
    if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
        for (int j = 1; j <= gv->NY; j++) {
            n = 1;
            for (int l = 1; l <= fdo; l++) {
                vy[j][gv->NX + l] = bufferlef_to_rig[j][n++];
            }
            for (int l = 1; l <= fdo - 1; l++) {
                vx[j][gv->NX + l] = bufferlef_to_rig[j][n++];
            }
        }

    /* no exchange if periodic boundary condition is applied */
    if ((gv->BOUNDARY) || (gv->POS[1] != 0))    /* no boundary exchange at left edge of global grid */
        for (int j = 1; j <= gv->NY; j++) {
            n = 1;
            for (int l = 1; l <= fdo - 1; l++) {
                vy[j][1 - l] = bufferrig_to_lef[j][n++];
            }
            for (int l = 1; l <= fdo; l++) {
                vx[j][1 - l] = bufferrig_to_lef[j][n++];
            }
        }
}
