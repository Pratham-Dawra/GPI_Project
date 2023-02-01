
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

/* ------------------------------------------------------------------------
 * prepare of update of the stress tensor
 * ------------------------------------------------------------------------*/

#include "fd.h"

void prepare_update_s(MemModel * mpm, GlobVar * gv)
{
    for (int l = 1; l <= gv->L; l++) {
        mpm->etajm[l] = mpm->peta[l];
        mpm->etaip[l] = mpm->peta[l];
    }
    for (int j = 1; j <= gv->NY; j++) {
        for (int i = 1; i <= gv->NX; i++) {
            mpm->fipjp[j][i] = mpm->puipjp[j][i] * gv->DT * (1.0 + gv->L * mpm->ptausipjp[j][i]);
            mpm->f[j][i] = mpm->pu[j][i] * gv->DT * (1.0 + gv->L * mpm->ptaus[j][i]);
            mpm->g[j][i] = mpm->ppi[j][i] * gv->DT * (1.0 + gv->L * mpm->ptaup[j][i]);
            for (int l = 1; l <= gv->L; l++) {
                mpm->bip[l] = 1.0 / (1.0 + (mpm->etaip[l] * 0.5));
                mpm->bjm[l] = 1.0 / (1.0 + (mpm->etajm[l] * 0.5));
                mpm->cip[l] = 1.0 - (mpm->etaip[l] * 0.5);
                mpm->cjm[l] = 1.0 - (mpm->etajm[l] * 0.5);
                mpm->dip[j][i][l] = mpm->puipjp[j][i] * mpm->etaip[l] * mpm->ptausipjp[j][i];
                mpm->d[j][i][l] = mpm->pu[j][i] * mpm->etajm[l] * mpm->ptaus[j][i];
                mpm->e[j][i][l] = mpm->ppi[j][i] * mpm->etajm[l] * mpm->ptaup[j][i];
            }
        }
    }
}
