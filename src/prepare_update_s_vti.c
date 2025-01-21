
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

void prepare_update_s_vti(MemModel * mpm, GlobVar * gv)
{
    for (int j = 1; j <= gv->NY; j++) {
        for (int i = 1; i <= gv->NX; i++) {
            /* unrelaxed moduli */
            mpm->pc55ipjpu[j][i] = mpm->pc55ipjp[j][i] * gv->DT * (1.0 + gv->L * mpm->ptau55ipjp[j][i]);
            mpm->pc13u[j][i] = mpm->pc13[j][i] * gv->DT * (1.0 + gv->L * mpm->ptau13[j][i]);
            mpm->pc11u[j][i] = mpm->pc11[j][i] * gv->DT * (1.0 + gv->L * mpm->ptau11[j][i]);
            mpm->pc33u[j][i] = mpm->pc33[j][i] * gv->DT * (1.0 + gv->L * mpm->ptau33[j][i]);
            for (int l = 1; l <= gv->L; l++) {
                mpm->bip[l] = 1.0 / (1.0 + (mpm->peta[l] * 0.5));
                mpm->cip[l] = 1.0 - (mpm->peta[l] * 0.5);
                /* module defects for each relaxation mechanism */
                mpm->pc55ipjpd[j][i][l] = mpm->pc55ipjp[j][i] * mpm->peta[l] * mpm->ptau55ipjp[j][i];
                mpm->pc13d[j][i][l] = mpm->pc13[j][i] * mpm->peta[l] * mpm->ptau13[j][i];
                mpm->pc33d[j][i][l] = mpm->pc33[j][i] * mpm->peta[l] * mpm->ptau33[j][i];
                mpm->pc11d[j][i][l] = mpm->pc11[j][i] * mpm->peta[l] * mpm->ptau11[j][i];
            }
        }
    }
}
