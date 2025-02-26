
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

/* ----------------------------------------------------------------------
 * Computation of local source coordinates
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include <math.h>
#include "logging.h"

float **splitsrc(float **srcpos, int *nsrc_loc, int nsrc, int *srcswitch, GlobVar *gv)
{
    int a, b, i = 0;
    float **srcpos_local = NULL;
    float **srcpos_dummy = matrix(1, NSPAR, 1, nsrc);

    for (int j = 1; j <= nsrc; j++) {
        a = (int)floor(srcpos[1][j] / gv->DH) + 1;  // convert x coordinate to index grid point
        b = (int)floor(srcpos[2][j] / gv->DH) + 1;  // convert y coordinate to index grid point
        log_info("j: %d, nsrc_loc: %d, nsrc: %d.\n", j, nsrc_loc, nsrc);
        srcswitch[j] = 0;
        if (a >= gv->GGRID[1] && a <= gv->GGRID[2] && b >= gv->GGRID[3] && b <= gv->GGRID[4]) {
            i++;
            srcpos_dummy[1][i] = (float)(a - gv->GGRID[1] + 1);
            srcpos_dummy[2][i] = (float)(b - gv->GGRID[3] + 1);
            srcpos_dummy[3][i] = 0.0f;
            srcpos_dummy[4][i] = srcpos[4][j];
            srcpos_dummy[5][i] = srcpos[5][j];
            srcpos_dummy[6][i] = srcpos[6][j];
            srcpos_dummy[7][i] = srcpos[7][j];
            srcpos_dummy[8][i] = srcpos[8][j];
            srcpos_dummy[9][i] = srcpos[9][j];
            srcpos_dummy[10][i] = srcpos[10][j];
            srcpos_dummy[11][i] = srcpos[11][j];
            srcpos_dummy[12][i] = srcpos[12][j];
            srcswitch[j] = 1;
        }
    }

    if (i > 0)
        srcpos_local = matrix(1, NSPAR, 1, i);
    for (int k = 1; k <= i; k++) {
        srcpos_local[1][k] = srcpos_dummy[1][k];
        srcpos_local[2][k] = srcpos_dummy[2][k];
        srcpos_local[3][k] = srcpos_dummy[3][k];
        srcpos_local[4][k] = srcpos_dummy[4][k];
        srcpos_local[5][k] = srcpos_dummy[5][k];
        srcpos_local[6][k] = srcpos_dummy[6][k];
        srcpos_local[7][k] = srcpos_dummy[7][k];
        srcpos_local[8][k] = srcpos_dummy[8][k];
        srcpos_local[9][k] = srcpos_dummy[9][k];
        srcpos_local[10][k] = srcpos_dummy[10][k];
        srcpos_local[11][k] = srcpos_dummy[11][k];
        srcpos_local[12][k] = srcpos_dummy[12][k];
    }
    free_matrix(srcpos_dummy, 1, NSPAR, 1, nsrc);

    *nsrc_loc = i;
    
    log_info("nsrc_loc: %d, nsrc: %d, srcswitch: %d.\n", *nsrc_loc, nsrc, srcswitch[1]);

    return srcpos_local;
}
