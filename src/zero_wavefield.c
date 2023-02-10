
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

void zero_wavefield(MemWavefield *mpw, GlobVar *gv)
{

    if (gv->L) {
        /* viscoelastic */
        for (int j = (-gv->ND + 1); j <= (gv->NY + gv->ND); j++) {
            for (int i = (-gv->ND + 1); i <= (gv->NX + gv->ND); i++) {
                zero_elastic(j, i, mpw);
                for (int l = 1; l <= gv->L; l++) {
                    zero_visco(j, i, l, mpw);
                }
            }
        }
    } else {
        /* elastic */
        for (int j = (-gv->ND + 1); j <= (gv->NY + gv->ND); j++) {
            for (int i = (-gv->ND + 1); i <= (gv->NX + gv->ND); i++) {
                zero_elastic(j, i, mpw);
            }
        }
    }

    if (gv->ABS_TYPE == 1) {
        /* PML Boundary */
        for (int j = 1; j <= gv->NY; j++) {
            for (int i = 1; i <= 2 * gv->FW; i++) {
                zero_PML_x(j, i, mpw);
            }
        }
        for (int j = 1; j <= 2 * gv->FW; j++) {
            for (int i = 1; i <= gv->NX; i++) {
                zero_PML_y(j, i, mpw);
            }
        }
    }
        
    if (gv->FDORDER_TIME == 4) {
        if (gv->L) {
        /* viscoelastic, FDORDER_TIME = 4 */
        for (int j = (-gv->ND + 1); j <= (gv->NY + gv->ND); j++) {
                for (int i = (-gv->ND + 1); i <= (gv->NX + gv->ND); i++) {
                    zero_elastic_4(j, i, mpw);
                    for (int l = 1; l <= gv->L; l++) {
                        zero_visco_4(j, i, l, mpw);
                    }
                }
            }
        } else {
            /* elastic, FDORDER_TIME = 4 */
            for (int j = (-gv->ND + 1); j <= (gv->NY + gv->ND); j++) {
                for (int i = (-gv->ND + 1); i <= (gv->NX + gv->ND); i++) {
                    zero_elastic_4(j, i, mpw);
                }
            }
        }
    }
}
