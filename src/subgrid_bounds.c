
/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/*---------------------------------------------------------------
calculate for-loop bounds for each process/subgrid

   [1] [2]                 [3] [4] GX
 [1]|---|-------------------|---|
    |   |                   |   |
 [2]|---|-------------------|---|
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
 [3]|---|-------------------|---|
    |   |                   |   |
 [4]|---|-------------------|---|

 GY

Width of Absorbing Boundary is gv->FW Gridpoints.
nx1,nx2,ny1,ny2 are first are the bounds for each subgrid without
absorbing boundaries.
In case of number of processes NPROCX=1 and NPROCY=1 all boundaries are in one
grid.
For NPROCX>2 and NPROCY>2 there are 9 different subgrids.
GX[2],GY[2] is the last gridpoint (FW) of the left boundary 
GX[3],GY[3] is the last gridpoint of the midpart  GX[3]=(nx2-FW),GY[3]=(ny2-FW)

Free Surface and Periodic Boundary conditions alter the bounds (see end of this file)
*--------------------------------------------------------------*/

#include "fd.h"

void subgrid_bounds(int nx1, int nx2, int ny1, int ny2, GlobVar * gv)
{
    /* GRID */
    switch (gv->NPROCY) {

      case 1:                  /*case gv->NPROCY=1 */

          switch (gv->NPROCX) {
            case 1:            /* Case gv->NPROCX=1 */
                gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;

                break;

            case 2:            /* Case gv->NPROCX=2 */
                if (gv->POS[1] == 0) {  /* left */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if (gv->POS[1] == gv->NPROCX - 1) { /* right */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                break;

            default:           /* Case gv->NPROCX>2 */

                if (gv->POS[1] == 0) {  /* left */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }

                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1)) {  /* mid */

                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if (gv->POS[1] == gv->NPROCX - 1) { /* right */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                break;
          }
          break;

      case 2:                  /*case gv->NPROCY=2 */
          switch (gv->NPROCX) {
            case 1:            /* Case gv->NPROCX=1 */
                if (gv->POS[2] == 0) {  /* top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if (gv->POS[2] == gv->NPROCY - 1) { /* bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                break;

            case 2:            /* Case gv->NPROCX=2 */
                if ((gv->POS[1] == 0) && (gv->POS[2] == 0)) {   /* corner left-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == 0) && (gv->POS[2] == gv->NPROCY - 1)) {  /* corner left-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == 0)) {  /* corner right-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) { /* corner right-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                break;

            default:           /* Case gv->NPROCX>2 */

                if ((gv->POS[1] == 0) && (gv->POS[2] == 0)) {   /* corner left-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == 0) && (gv->POS[2] == gv->NPROCY - 1)) {  /* corner left-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == 0)) {  /* corner right-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) { /* corner right-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1) && (gv->POS[2] == 0)) { /*top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) {    /* bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }

                break;
          }                     /* end of switch gv->NPROCX */
          break;

      default:                 /*case gv->NPROCY>2 */
          switch (gv->NPROCX) {
            case 1:            /* Case gv->NPROCX=1 */
                if (gv->POS[2] == 0) {  /* top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                if ((gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) {  /* mid */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                if (gv->POS[2] == gv->NPROCY - 1) { /* bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                break;

            case 2:            /* Case gv->NPROCX=2 */
                if ((gv->POS[1] == 0) && (gv->POS[2] == 0)) {   /* corner left-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == 0) && (gv->POS[2] == gv->NPROCY - 1)) {  /* corner left-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == 0)) {  /* corner right-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) { /* corner right-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if ((gv->POS[1] == 0) && (gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) { /* left */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) {    /* right */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                break;

            default:           /* Case gv->NPROCX>2 */

                if ((gv->POS[1] == 0) && (gv->POS[2] == 0)) {   /* corner left-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == 0) && (gv->POS[2] == gv->NPROCY - 1)) {  /* corner left-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == 0)) {  /* corner right-top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) { /* corner right-bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                if ((gv->POS[1] == 0) && (gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) { /* left */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = gv->FW, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] == gv->NPROCX - 1) && (gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) {    /* right */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2 - gv->FW, gv->GX[4] = nx2;
                }

                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1) && (gv->POS[2] == 0)) { /*top */
                    gv->GY[1] = ny1, gv->GY[2] = gv->FW, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1) && (gv->POS[2] == gv->NPROCY - 1)) {    /* bottom */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2 - gv->FW, gv->GY[4] = ny2;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }

                if ((gv->POS[1] != 0) && (gv->POS[1] != gv->NPROCX - 1) && (gv->POS[2] != 0) && (gv->POS[2] != gv->NPROCY - 1)) {   /* center */
                    gv->GY[1] = ny1, gv->GY[2] = ny1 - 1, gv->GY[3] = ny2, gv->GY[4] = ny2 - 1;
                    gv->GX[1] = nx1, gv->GX[2] = nx1 - 1, gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
                }
                break;
          }                     /*end of switch (gv->NPROCX) */
          break;
    }                           /*end of switch (gv->NPROCY) */

    if (gv->FREE_SURF == 1)
        if (gv->POS[2] == 0)
            gv->GY[2] = ny1 - 1;

    /* Periodic Boundary condition */
    if (gv->BOUNDARY == 1) {
        if (gv->POS[1] == 0)
            gv->GX[2] = nx1 - 1;
        if (gv->POS[1] == gv->NPROCX - 1)
            gv->GX[3] = nx2, gv->GX[4] = nx2 - 1;
    }
}
