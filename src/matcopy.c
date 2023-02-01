
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
 * For the averaging of material properties each process requires values
 * at the indices 0 and NX+1 etc. These lie on the neighbouring processes.
 * Thus, they have to be copied which is done by this function.
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void matcopy(float **rho, float **pi, float **u, float **taus, float **taup, GlobVar * gv)
{
    MPI_Status status;
    double time1, time2;
    float **bufferlef_to_rig, **bufferrig_to_lef;
    float **buffertop_to_bot, **bufferbot_to_top;

    bufferlef_to_rig = matrix(0, gv->NY + 1, 1, 5);
    bufferrig_to_lef = matrix(0, gv->NY + 1, 1, 5);
    buffertop_to_bot = matrix(0, gv->NX + 1, 1, 5);
    bufferbot_to_top = matrix(0, gv->NX + 1, 1, 5);

    log_debugc(0, "Copying material properties at inner boundaries.\n");
    time1 = MPI_Wtime();

    //    if (gv->POS[2]!=0)  /* no boundary exchange at top of global grid */
    for (int i = 0; i <= gv->NX + 1; i++) {
        /* storage of top of local volume into buffer */
        buffertop_to_bot[i][1] = rho[1][i];
        buffertop_to_bot[i][2] = pi[1][i];
        buffertop_to_bot[i][3] = u[1][i];
        buffertop_to_bot[i][4] = taus[1][i];
        buffertop_to_bot[i][5] = taup[1][i];
    }

    //    if (gv->POS[2]!=NPROCY-1)   /* no boundary exchange at bottom of global grid */
    for (int i = 0; i <= gv->NX + 1; i++) {
        /* storage of bottom of local volume into buffer */
        bufferbot_to_top[i][1] = rho[gv->NY][i];
        bufferbot_to_top[i][2] = pi[gv->NY][i];
        bufferbot_to_top[i][3] = u[gv->NY][i];
        bufferbot_to_top[i][4] = taus[gv->NY][i];
        bufferbot_to_top[i][5] = taup[gv->NY][i];
    }

  /*=========sending and receiving of the boundaries==========*/

    MPI_Bsend(&buffertop_to_bot[1][1], gv->NX * 5, MPI_FLOAT, gv->INDEX[3], gv->TAG5, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Recv(&buffertop_to_bot[1][1], gv->NX * 5, MPI_FLOAT, gv->INDEX[4], gv->TAG5, MPI_COMM_WORLD, &status);
    MPI_Bsend(&bufferbot_to_top[1][1], gv->NX * 5, MPI_FLOAT, gv->INDEX[4], gv->TAG6, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Recv(&bufferbot_to_top[1][1], gv->NX * 5, MPI_FLOAT, gv->INDEX[3], gv->TAG6, MPI_COMM_WORLD, &status);

    //    if (gv->POS[2]!=NPROCY-1)   /* no boundary exchange at bottom of global grid */
    for (int i = 0; i <= gv->NX + 1; i++) {
        rho[gv->NY + 1][i] = buffertop_to_bot[i][1];
        pi[gv->NY + 1][i] = buffertop_to_bot[i][2];
        u[gv->NY + 1][i] = buffertop_to_bot[i][3];
        taus[gv->NY + 1][i] = buffertop_to_bot[i][4];
        taup[gv->NY + 1][i] = buffertop_to_bot[i][5];
    }

    //    if (gv->POS[2]!=0)  /* no boundary exchange at top of global grid */
    for (int i = 0; i <= gv->NX + 1; i++) {
        rho[0][i] = bufferbot_to_top[i][1];
        pi[0][i] = bufferbot_to_top[i][2];
        u[0][i] = bufferbot_to_top[i][3];
        taus[0][i] = bufferbot_to_top[i][4];
        taup[0][i] = bufferbot_to_top[i][5];
    }

    //    if (gv->POS[1]!=0)  /* no boundary exchange at left edge of global grid */
    for (int j = 0; j <= gv->NY + 1; j++) {
        /* storage of left edge of local volume into buffer */
        bufferlef_to_rig[j][1] = rho[j][1];
        bufferlef_to_rig[j][2] = pi[j][1];
        bufferlef_to_rig[j][3] = u[j][1];
        bufferlef_to_rig[j][4] = taus[j][1];
        bufferlef_to_rig[j][5] = taup[j][1];
    }

    //    if (gv->POS[1]!=NPROCX-1)   /* no boundary exchange at right edge of global grid */
    for (int j = 0; j <= gv->NY + 1; j++) {
        /* storage of right edge of local volume into buffer */
        bufferrig_to_lef[j][1] = rho[j][gv->NX];
        bufferrig_to_lef[j][2] = pi[j][gv->NX];
        bufferrig_to_lef[j][3] = u[j][gv->NX];
        bufferrig_to_lef[j][4] = taus[j][gv->NX];
        bufferrig_to_lef[j][5] = taup[j][gv->NX];
    }

    MPI_Bsend(&bufferlef_to_rig[0][1], (gv->NY + 2) * 5, MPI_FLOAT, gv->INDEX[1], gv->TAG1, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Recv(&bufferlef_to_rig[0][1], (gv->NY + 2) * 5, MPI_FLOAT, gv->INDEX[2], gv->TAG1, MPI_COMM_WORLD, &status);
    MPI_Bsend(&bufferrig_to_lef[0][1], (gv->NY + 2) * 5, MPI_FLOAT, gv->INDEX[2], gv->TAG2, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Recv(&bufferrig_to_lef[0][1], (gv->NY + 2) * 5, MPI_FLOAT, gv->INDEX[1], gv->TAG2, MPI_COMM_WORLD, &status);

    //    if (gv->POS[1]!=NPROCX-1)   /* no boundary exchange at right edge of global grid */
    for (int j = 0; j <= gv->NY + 1; j++) {
        rho[j][gv->NX + 1] = bufferlef_to_rig[j][1];
        pi[j][gv->NX + 1] = bufferlef_to_rig[j][2];
        u[j][gv->NX + 1] = bufferlef_to_rig[j][3];
        taus[j][gv->NX + 1] = bufferlef_to_rig[j][4];
        taup[j][gv->NX + 1] = bufferlef_to_rig[j][5];
    }

    //    if (gv->POS[1]!=0)  /* no boundary exchange at left edge of global grid */
    for (int j = 0; j <= gv->NY + 1; j++) {
        rho[j][0] = bufferrig_to_lef[j][1];
        pi[j][0] = bufferrig_to_lef[j][2];
        u[j][0] = bufferrig_to_lef[j][3];
        taus[j][0] = bufferrig_to_lef[j][4];
        taup[j][0] = bufferrig_to_lef[j][5];
    }

    time2 = MPI_Wtime();
    log_debugc(0, "Finished copying material properties at inner boundaries (real time: %4.2fs).\n", time2 - time1);

    free_matrix(bufferlef_to_rig, 0, gv->NY + 1, 1, 5);
    free_matrix(bufferrig_to_lef, 0, gv->NY + 1, 1, 5);
    free_matrix(buffertop_to_bot, 0, gv->NX + 1, 1, 5);
    free_matrix(bufferbot_to_top, 0, gv->NX + 1, 1, 5);
}
