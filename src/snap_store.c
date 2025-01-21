
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

/*----------------------------------------------------------------------
 * Store 2D snapshot of current timestep for gradient calculation
 *----------------------------------------------------------------------*/

#include "fd.h"

void snap_store(int nt, int hin, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    int i, j;

    for (j = 1; j <= gv->NY; j++) {
        for (i = 1; i <= gv->NX; i++) {
            minv->forward_prop_rho_x[j][i][hin] = minv->pvxp1[j][i];
            minv->forward_prop_rho_y[j][i][hin] = minv->pvyp1[j][i];
            if (vinv->VELOCITY == 0) {
                minv->forward_prop_x[j][i][hin] = mpw->psxx[j][i];
                minv->forward_prop_y[j][i][hin] = mpw->psyy[j][i];
                minv->forward_prop_u[j][i][hin] = mpw->psxy[j][i];
            } else {
                minv->forward_prop_x[j][i][hin] = minv->ux[j][i];
                minv->forward_prop_y[j][i][hin] = minv->uy[j][i];
                minv->forward_prop_u[j][i][hin] = minv->uxy[j][i];
            }
        }
    }
    minv->DTINV_help[nt] = 1;
}
