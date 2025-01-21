
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
 *  Precalculate derivates of particle velocities
 * -------------------------------------------------------------*/

#include "fd.h"

void v_derivatives(MemWavefield * mpw, GlobVar * gv)
{
    float vxx, vyy, vxy, vyx;

    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {  /* elastic cases */
        for (int j = 1; j <= gv->NY; j++) {
            for (int i = 1; i <= gv->NX; i++) {
                gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
                mpw->pvxx[j][i] = vxx;
                mpw->pvyy[j][i] = vyy;
                mpw->pvyx[j][i] = vyx;
                mpw->pvxy[j][i] = vxy;
            }
        }
    } else {    /* acoustic cases */
        for (int j = 1; j <= gv->NY; j++) {
            for (int i = 1; i <= gv->NX; i++) {
                gv->FDOP_AC_S(i, j, &vxx, &vyy, mpw);
                mpw->pvxx[j][i] = vxx;
                mpw->pvyy[j][i] = vyy;
            }
        }
    }
}
