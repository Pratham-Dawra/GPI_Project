
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   zero wavefield
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_PML_elastic(int nx1, int nx2, int ny1, int ny2, MemWavefield * mpw, GlobVar * gv)
{
    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            mpw->pvx[j][i] = 0.0f;
            mpw->pvy[j][i] = 0.0f;
            mpw->psxx[j][i] = 0.0f;
            mpw->psyy[j][i] = 0.0f;
            mpw->psxy[j][i] = 0.0f;
        }
    }

    for (int j = 1; j <= gv->NY; j++) {
        for (int i = 1; i <= 2 * gv->FW; i++) {
            mpw->psi_sxx_x[j][i] = 0.0f;
            mpw->psi_sxy_x[j][i] = 0.0f;
            mpw->psi_vxx[j][i] = 0.0f;
            mpw->psi_vyx[j][i] = 0.0f;
            mpw->psi_vxxs[j][i] = 0.0f;
        }
    }

    for (int j = 1; j <= 2 * gv->FW; j++) {
        for (int i = 1; i <= gv->NX; i++) {
            mpw->psi_syy_y[j][i] = 0.0f;
            mpw->psi_sxy_y[j][i] = 0.0f;
            mpw->psi_vyy[j][i] = 0.0f;
            mpw->psi_vxy[j][i] = 0.0f;
        }
    }
}
