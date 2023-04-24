
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

/*
 * Update Function of the particle velocity wavefields 
 */

#include "fd.h"

void wavefield_update_v(int i, int j, int sw, float sxx_x, float sxy_x, float sxy_y, float syy_y, MemModel *mpm,
                        MemWavefield *mpw, MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{
    float dtdh = gv->DT / gv->DH;

    if (gv->MODE == FWI) {
        if (sw == 1) {  /* Backpropagation (sw==1) */
            if (vinv->VELOCITY == 0) {
                minv->pvxp1[j][i] += mpw->pvx[j][i] * gv->DT;
                minv->pvyp1[j][i] += mpw->pvy[j][i] * gv->DT;
            } else {
                minv->pvxp1[j][i] = mpw->pvx[j][i];
                minv->pvyp1[j][i] = mpw->pvy[j][i];
            }
        } else {    /* Forward Modelling (sw==0) */
            minv->pvxp1[j][i] = mpm->prip[j][i] * (sxx_x + sxy_y) / gv->DH;
            minv->pvyp1[j][i] = mpm->prjp[j][i] * (sxy_x + syy_y) / gv->DH;
        }
    }

    mpw->pvx[j][i] += (sxx_x + sxy_y) * dtdh * mpm->prip[j][i];
    mpw->pvy[j][i] += (sxy_x + syy_y) * dtdh * mpm->prjp[j][i];

}
