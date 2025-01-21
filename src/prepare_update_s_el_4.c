
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

/* ------------------------------------------------------------------------
 * ATTENTION: The parameters below will be scaled by factor c* due to
 * Adams-Bashforth method, so be aware and only call this function when
 * FDORDER_TIME is set to 4
 * ------------------------------------------------------------------------*/

#include "fd.h"

void prepare_update_s_el_4(MemModel *mpm, GlobVar *gv)
{

    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;

    for (int j = 1; j <= gv->NY; j++) {
        for (int i = 1; i <= gv->NX; i++) {
            float fipjp = mpm->puipjp[j][i] * gv->DT;
            float f = mpm->pu[j][i] * gv->DT;
            float g = mpm->ppi[j][i] * gv->DT;
        }
    }
}
