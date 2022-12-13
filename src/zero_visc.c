
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

/*------------------------------------------------------------------------
 * Initialization of the wave field with zero values (zero wavefield)
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_visc(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx, float **syy, float **sxy,
               float ***pr, float ***pp, float ***pq, GlobVar *gv)
{
    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            vx[j][i] = 0.0f;
            vy[j][i] = 0.0f;
            sxx[j][i] = 0.0f;
            syy[j][i] = 0.0f;
            sxy[j][i] = 0.0f;
        }
    }

    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            for (int l = 1; l <= gv->L; l++) {
                pr[j][i][l] = 0.0f;
                pp[j][i][l] = 0.0f;
                pq[j][i][l] = 0.0f;
            }
        }
    }
}
