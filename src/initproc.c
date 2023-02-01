
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

void initproc(GlobVar * gv)
{
    if ((gv->NPROC != gv->NP) && (gv->MPID == 0)) {
        log_error("You specified NPROCX*NPROCY=%d (in json file) but number of MPI processes set to NP=%d.\n",
                  gv->NPROC, gv->NP);
        log_fatal("Number of MPI processes (NP) and number of processing elements (NPROCX*NPROCY) differ!\n");
    }

    /* determine the length of the subarray on this processor */
    gv->IENDX = gv->NX / gv->NPROCX;
    gv->IENDY = gv->NY / gv->NPROCY;

    /* POS(1) indicates x POSition of the processor in the logical 3D processor array */
    if ((gv->NX % gv->NPROCX) > 0)
        log_fatal("NX%NPROCX (modulus) must be zero!\n");
    if ((gv->NY % gv->NPROCY) > 0)
        log_fatal("NY%NPROCY (modulus) must be zero!\n");

    log_debugc(0, "Size of subdomains (in grid points): %d (horz.) times %d (vert.)\n", gv->IENDX, gv->IENDY);

    /* index is indicating neighbouring processes */
    gv->INDEX[1] = gv->MPID - 1;    /* left    */
    gv->INDEX[2] = gv->MPID + 1;    /* right   */
    gv->INDEX[3] = gv->MPID - gv->NPROCX;   /* upper   */
    gv->INDEX[4] = gv->MPID + gv->NPROCX;   /* lower   */

    /* POS indicates the processor location in the 3D logical processor array */
    gv->POS[1] = gv->MPID % gv->NPROCX; /* x coordinate */
    gv->POS[2] = (gv->MPID / gv->NPROCX);   /* y coordinate */

    if (gv->POS[1] == 0)
        gv->INDEX[1] = gv->INDEX[1] + gv->NPROCX;
    if (gv->POS[1] == gv->NPROCX - 1)
        gv->INDEX[2] = gv->INDEX[2] - gv->NPROCX;
    if (gv->POS[2] == 0)
        gv->INDEX[3] = (gv->NPROCX * gv->NPROCY) + gv->MPID - gv->NPROCX;
    if (gv->POS[2] == gv->NPROCY - 1)
        gv->INDEX[4] = gv->MPID + gv->NPROCX - (gv->NPROCX * gv->NPROCY);

    log_debug("PE %d; col %d; row %d; left-right neighbor %d,%d; upper-lower neighbor %d,%d\n",
              gv->MPID, gv->POS[1], gv->POS[2], gv->INDEX[1], gv->INDEX[2], gv->INDEX[3], gv->INDEX[4]);

    return;
}

/* Example: NPROX=5, NPROCY=4 => NPROC=20           */

/*                                                  */

/* ================================================ */

/*                                                  */

/*  0               0: 4,1                  0: 15,5 */

/*  1               1: 0,2                  0: 16,6 */

/*  2               2: 1,3                  0: 17,7 */

/*  3               3: 2,4                  0: 18,8 */

/*  4               4: 3,0                  0: 19,9 */

/*  5               0: 9,6                  1: 0,10 */

/*  6               1: 5,7                  1: 1,11 */

/*  7               2: 6,8                  1: 2,12 */

/*  8               3: 7,9                  1: 3,13 */

/*  9               4: 8,5                  1: 4,14 */

/*  10              0: 14,11                2: 5,15 */

/*  11              1: 10,12                2: 6,16 */

/*  12              2: 11,13                2: 7,17 */

/*  13              3: 12,14                2: 8,18 */

/*  14              4: 13,10                2: 9,19 */

/*  15              0: 19,16                3: 10,0 */

/*  16              1: 15,17                3: 11,1 */

/*  17              2: 16,18                3: 12,2 */

/*  18              3: 17,19                3: 13,3 */

/*  19              4: 18,15                3: 14,4 */

/*                                                  */

/* ================================================ */

/*                                                  */

/*       0         1        2       3       4       */

/*       x         x        x       x       x       */

/*                                                  */

/*       5         6        7       8       9       */

/*       x         x        x       x       x       */

/*                                                  */

/*      10        11       12      13      14       */

/*       x         x        x       x       x       */

/*                                                  */

/*      15        16       17      18      19       */

/*       x         x        x       x       x       */

/*                                                  */

/* ================================================ */
