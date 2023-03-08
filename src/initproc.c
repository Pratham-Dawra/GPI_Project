
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
 * This is function initproc.
 * Dividing the 2-D FD grid into domains and assigning the
 * PEs to these domains,
 * -------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void initproc(GlobVar *gv)
{
    int nmpi;
    MPI_Comm_size(MPI_COMM_WORLD, &nmpi);
    if ((gv->NPROC != nmpi) && (gv->MPID == 0)) {
        log_error("You specified NPROCX*NPROCY=%d (in json file) but number of MPI processes set to NP=%d.\n",
                  gv->NPROC, nmpi);
        log_fatal("Number of MPI processes (NP) and number of processing elements (NPROCX*NPROCY) differ!\n");
    }

    /* determine the length of the subarray on this processor */
    gv->NX = gv->NXG / gv->NPROCX;
    gv->NY = gv->NYG / gv->NPROCY;

    if ((gv->NXG % gv->NPROCX) > 0)
        log_fatal("NX%NPROCX (modulus) must be zero!\n");
    if ((gv->NYG % gv->NPROCY) > 0)
        log_fatal("NY%NPROCY (modulus) must be zero!\n");

    log_debugc(0, "Size of subdomains (in grid points): %d (horz.) times %d (vert.)\n", gv->NX, gv->NY);

    /* index is indicating neighbouring processes */
    gv->INDEX[1] = gv->MPID - 1;    /* left    */
    gv->INDEX[2] = gv->MPID + 1;    /* right   */
    gv->INDEX[3] = gv->MPID - gv->NPROCX;   /* upper   */
    gv->INDEX[4] = gv->MPID + gv->NPROCX;   /* lower   */

    /* POS indicates the processor location in the 3D logical processor array (x=1, y=2) */
    gv->POS[1] = gv->MPID % gv->NPROCX;     /* x coordinate */
    gv->POS[2] = (gv->MPID / gv->NPROCX);   /* y coordinate */

    if (gv->POS[1] == 0)
        gv->INDEX[1] = gv->INDEX[1] + gv->NPROCX;
    if (gv->POS[1] == gv->NPROCX - 1)
        gv->INDEX[2] = gv->INDEX[2] - gv->NPROCX;
    if (gv->POS[2] == 0)
        gv->INDEX[3] = (gv->NPROCX * gv->NPROCY) + gv->MPID - gv->NPROCX;
    if (gv->POS[2] == gv->NPROCY - 1)
        gv->INDEX[4] = gv->MPID + gv->NPROCX - (gv->NPROCX * gv->NPROCY);

    /* GGRID are the global model indices [1,NYG][1,NXG] handled by this logical processor, *
     * 1=left (idx xmin), 2=right (idx xmax), 3=upper (idx ymin), 4=lower (idx ymax)        */
    gv->GGRID[1] = 1;
    for (int i = 0; i < gv->POS[1]; ++i) {
        gv->GGRID[1] += gv->NX;
    }
    gv->GGRID[2] = gv->GGRID[1] + (gv->NX - 1);
    gv->GGRID[3] = 1;
    for (int i = 0; i < gv->POS[2]; ++i) {
        gv->GGRID[3] += gv->NY;
    }
    gv->GGRID[4] = gv->GGRID[3] + (gv->NY - 1);

    log_debug("PE %d; col %d; row %d; neighbors l-r %d,%d, u-l %d,%d; grid points xb-xe %d,%d yb-ye %d,%d\n",
              gv->MPID, gv->POS[1], gv->POS[2], gv->INDEX[1], gv->INDEX[2], gv->INDEX[3], gv->INDEX[4],
              gv->GGRID[1], gv->GGRID[2], gv->GGRID[3], gv->GGRID[4]);

    /* if IDX and/or IDY not equal to 1, snapshots are coarser than the model grid; find *
     * global indices each logical processor has to output;  1=left (idx xmin), ...      */
    if (gv->IDX == 1) {
        gv->SNAPIDX[1] = gv->GGRID[1];
        gv->SNAPIDX[2] = gv->GGRID[2];
    } else {
        gv->SNAPIDX[1] = gv->NXG+1;
        gv->SNAPIDX[2] = -1;
        for (int i = 1; i <= gv->NXG; i += gv->IDX) {
            if (i >= gv->GGRID[1] && i <= gv->GGRID[2]) {
                if (i < gv->SNAPIDX[1]) gv->SNAPIDX[1] = i;
                if (i > gv->SNAPIDX[2]) gv->SNAPIDX[2] = i;
            }
        }
    }
    if (gv->IDY == 1) {
        gv->SNAPIDX[3] = gv->GGRID[3];
        gv->SNAPIDX[4] = gv->GGRID[4];
    } else {
        gv->SNAPIDX[3] = gv->NYG+1;
        gv->SNAPIDX[4] = -1;
        for (int j = 1; j <= gv->NYG; j += gv->IDY) {
            if (j >= gv->GGRID[3] && j <= gv->GGRID[4]) {
                if (j < gv->SNAPIDX[3]) gv->SNAPIDX[3] = j;
                if (j > gv->SNAPIDX[4]) gv->SNAPIDX[4] = j;
            }
        }
    }
    /* now convert global indices to local processor indices */
    gv->SNAPIDX[1] = gv->SNAPIDX[1] - gv->GGRID[1] + 1;
    gv->SNAPIDX[2] = gv->SNAPIDX[2] - gv->GGRID[1] + 1;
    gv->SNAPIDX[3] = gv->SNAPIDX[3] - gv->GGRID[3] + 1;
    gv->SNAPIDX[4] = gv->SNAPIDX[4] - gv->GGRID[3] + 1;

    log_debug("PE %d; col %d; row %d; snap indices xs-xe %d,%d yb-ye %d,%d\n",
              gv->MPID, gv->POS[1], gv->POS[2], gv->SNAPIDX[1], gv->SNAPIDX[2], gv->SNAPIDX[3], gv->SNAPIDX[4]);

    return;
}

/* Example: NPROCX=5, NPROCY=4 => NPROC=20          *
 * ================================================ *
 * Rank=MPID, x=POS[1], y=POS[2]                    *
 * l=INDEX[1], r=INDEX[2], u=INDEX[3], l=INDEX[4]   *
 * ================================================ *
 *  Rank            x: l,r             y: u,l       *
 *  ----            -------            -------      *
 *  0               0: 4,1             0: 15,5      *
 *  1               1: 0,2             0: 16,6      *
 *  2               2: 1,3             0: 17,7      *
 *  3               3: 2,4             0: 18,8      *
 *  4               4: 3,0             0: 19,9      *
 *  5               0: 9,6             1: 0,10      *
 *  6               1: 5,7             1: 1,11      *
 *  7               2: 6,8             1: 2,12      *
 *  8               3: 7,9             1: 3,13      *
 *  9               4: 8,5             1: 4,14      *
 *  10              0: 14,11           2: 5,15      *
 *  11              1: 10,12           2: 6,16      *
 *  12              2: 11,13           2: 7,17      *
 *  13              3: 12,14           2: 8,18      *
 *  14              4: 13,10           2: 9,19      *
 *  15              0: 19,16           3: 10,0      *
 *  16              1: 15,17           3: 11,1      *
 *  17              2: 16,18           3: 12,2      *
 *  18              3: 17,19           3: 13,3      *
 *  19              4: 18,15           3: 14,4      *
 *                                                  *
 * ================================================ *
 *                                                  *
 *       0         1        2       3       4       *
 *       x         x        x       x       x       *
 *                                                  *
 *       5         6        7       8       9       *
 *       x         x        x       x       x       *
 *                                                  *
 *      10        11       12      13      14       *
 *       x         x        x       x       x       *
 *                                                  *
 *      15        16       17      18      19       *
 *       x         x        x       x       x       *
 *                                                  *
 * ================================================ */
