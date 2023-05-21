
/*------------------------------------------------------------------------
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
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 * Damping functions for the stress fields and the particle velocity fields
 * used at the absorbing boundaries (ABS=2).    
 * ----------------------------------------------------------------------*/

#include "fd.h"

void abs_update_s_ac2(int i, int j, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psxx[j][i] *= mpm->absorb_coeff[j][i];
}

void abs_update_s_ac1(int i, int j, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psxx[j][i] *= mpm->absorb_coeff[j][i];
    mpw->psyy[j][i] *= mpm->absorb_coeff[j][i];

}

void abs_update_s(int i, int j, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psxy[j][i] *= mpm->absorb_coeff[j][i];
    mpw->psyy[j][i] *= mpm->absorb_coeff[j][i];
    mpw->psxx[j][i] *= mpm->absorb_coeff[j][i];
}

void abs_update_v(int i, int j, MemModel * mpm, MemWavefield * mpw)
{
    mpw->pvx[j][i] *= mpm->absorb_coeff[j][i];
    mpw->pvy[j][i] *= mpm->absorb_coeff[j][i];
}
