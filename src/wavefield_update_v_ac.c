
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

void wavefield_update_v_ac(int i, int j, float sxx_x, float sxx_y, MemModel * mpm,
                        MemWavefield * mpw, GlobVar * gv)
{
    float dtdh = gv->DT / gv->DH;
    mpw->pvx[j][i] += (sxx_x) * dtdh * mpm->prip[j][i];
    mpw->pvy[j][i] += (syy_y) * dtdh * mpm->prjp[j][i];
}
