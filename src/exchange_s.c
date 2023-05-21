
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
 *   local grid into buffer arrays and  exchanged between processes.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_s(MemWavefield * mpw, GlobVar * gv)
{
    MPI_Status status;
    int n;

    int fdo = gv->ND;
    int fdo2 = 2 * fdo;

    
    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {  /* elastic cases */
        
        /* top - bottom */
        
        if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                /* storage of top of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->buffertop_to_bot[i][n++] = mpw->psxy[l][i];
                }
                for (int l = 1; l <= fdo; l++) {
                    mpw->buffertop_to_bot[i][n++] = mpw->psyy[l][i];
                }
            }
        
        if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                /* storage of bottom of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->bufferbot_to_top[i][n++] = mpw->psxy[gv->NY - l + 1][i];
                }
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->bufferbot_to_top[i][n++] = mpw->psyy[gv->NY - l + 1][i];
                }
            }
        /* send and receive values for points at inner boundaries */
        MPI_Sendrecv_replace(&mpw->buffertop_to_bot[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[3], gv->TAG5, gv->INDEX[4],
                             gv->TAG5, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(&mpw->bufferbot_to_top[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[4], gv->TAG6, gv->INDEX[3],
                             gv->TAG6, MPI_COMM_WORLD, &status);
        
        if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psxy[gv->NY + l][i] = mpw->buffertop_to_bot[i][n++];
                }
                for (int l = 1; l <= fdo; l++) {
                    mpw->psyy[gv->NY + l][i] = mpw->buffertop_to_bot[i][n++];
                }
            }
        
        if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->psxy[1 - l][i] = mpw->bufferbot_to_top[i][n++];
                }
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psyy[1 - l][i] = mpw->bufferbot_to_top[i][n++];
                }
            }
        
        /* left - right */
        
        if ((gv->BOUNDARY) || (gv->POS[1] != 0))    /* no boundary exchange at left edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                /* storage of left edge of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->bufferlef_to_rig[j][n++] = mpw->psxy[j][l];
                }
                for (int l = 1; l <= fdo; l++) {
                    mpw->bufferlef_to_rig[j][n++] = mpw->psxx[j][l];
                }
            }
        
        if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                /* storage of right edge of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->bufferrig_to_lef[j][n++] = mpw->psxy[j][gv->NX - l + 1];
                }
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->bufferrig_to_lef[j][n++] = mpw->psxx[j][gv->NX - l + 1];
                }
            }
        
        /* send and receive values for points at inner boundaries */
        MPI_Sendrecv_replace(&mpw->bufferlef_to_rig[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[1], gv->TAG1, gv->INDEX[2],
                             gv->TAG1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(&mpw->bufferrig_to_lef[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[2], gv->TAG2, gv->INDEX[1],
                             gv->TAG2, MPI_COMM_WORLD, &status);
        /* send and reveive values at edges of the local grid */
        
        if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psxy[j][gv->NX + l] = mpw->bufferlef_to_rig[j][n++];
                }
                for (int l = 1; l <= fdo; l++) {
                    mpw->psxx[j][gv->NX + l] = mpw->bufferlef_to_rig[j][n++];
                }
            }
        
        if ((gv->BOUNDARY) || (gv->POS[1] != 0))    /* no boundary exchange at left edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->psxy[j][1 - l] = mpw->bufferrig_to_lef[j][n++];
                }
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psxx[j][1 - l] = mpw->bufferrig_to_lef[j][n++];
                }
            }
        
    }
    else /* acoustic cases */
    {
        
        /*fdo2 = fdo;*/
        /* top - bottom */
        
        if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                /* storage of top of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->buffertop_to_bot[i][n++] = mpw->psyy[l][i];
                }
            }
        
        if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                /* storage of bottom of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->bufferbot_to_top[i][n++] = mpw->psyy[gv->NY - l + 1][i];
                }
            }
        /* send and receive values for points at inner boundaries */
        MPI_Sendrecv_replace(&mpw->buffertop_to_bot[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[3], gv->TAG5, gv->INDEX[4],
                             gv->TAG5, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(&mpw->bufferbot_to_top[1][1], gv->NX * fdo2, MPI_FLOAT, gv->INDEX[4], gv->TAG6, gv->INDEX[3],
                             gv->TAG6, MPI_COMM_WORLD, &status);
        
        if (gv->POS[2] != gv->NPROCY - 1)   /* no boundary exchange at bottom of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                n = 1;
                 for (int l = 1; l <= fdo; l++) {
                    mpw->psyy[gv->NY + l][i] = mpw->buffertop_to_bot[i][n++];
                }
            }
        
        if (gv->POS[2] != 0)        /* no boundary exchange at top of global grid */
            for (int i = 1; i <= gv->NX; i++) {
                n = 1;
                 for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psyy[1 - l][i] = mpw->bufferbot_to_top[i][n++];
                }
            }
        
        /* left - right */
        
        if ((gv->BOUNDARY) || (gv->POS[1] != 0))    /* no boundary exchange at left edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                /* storage of left edge of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo; l++) {
                    mpw->bufferlef_to_rig[j][n++] = mpw->psxx[j][l];
                }
            }
        
        if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                /* storage of right edge of local volume into buffer */
                n = 1;
                for (int l = 1; l <= fdo - 1; l++) {
                    mpw->bufferrig_to_lef[j][n++] = mpw->psxx[j][gv->NX - l + 1];
                }
            }
        
        /* send and receive values for points at inner boundaries */
        MPI_Sendrecv_replace(&mpw->bufferlef_to_rig[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[1], gv->TAG1, gv->INDEX[2],
                             gv->TAG1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(&mpw->bufferrig_to_lef[1][1], gv->NY * fdo2, MPI_FLOAT, gv->INDEX[2], gv->TAG2, gv->INDEX[1],
                             gv->TAG2, MPI_COMM_WORLD, &status);
        /* send and reveive values at edges of the local grid */
        
        if ((gv->BOUNDARY) || (gv->POS[1] != gv->NPROCX - 1))   /* no boundary exchange at right edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                n = 1;
                 for (int l = 1; l <= fdo; l++) {
                    mpw->psxx[j][gv->NX + l] = mpw->bufferlef_to_rig[j][n++];
                }
            }
        
        if ((gv->BOUNDARY) || (gv->POS[1] != 0))    /* no boundary exchange at left edge of global grid */
            for (int j = 1; j <= gv->NY; j++) {
                n = 1;
                 for (int l = 1; l <= fdo - 1; l++) {
                    mpw->psxx[j][1 - l] = mpw->bufferrig_to_lef[j][n++];
                }
            }
        
        
    }
}
