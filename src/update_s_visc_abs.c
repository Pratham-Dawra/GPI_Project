
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
 *   updating stress components at gridpoints of the absorbing frame (ABS=2 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_abs(int *gx, int *gy,
                       float **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float ***r, float ***p, float ***q,
                       float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                       float *cjm, float ***d, float ***e, float ***dip, float **absorb_coeff, GlobVar *gv)
{
    float vxx, vyy, vxy, vyx;

    for (int j = gy[2] + 1; j <= gy[3]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* right boundary */
    for (int j = gy[2] + 1; j <= gy[3]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* top boundary */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[2] + 1; i <= gx[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* bottom boundary */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[2] + 1; i <= gx[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* corners */

    /*left-top */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /*left-bottom */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* right-top */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }

    /* right-bottom */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            wavefield_update_s_visc(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, r, p, q,
                                    fipjp, f, g, bip, bjm, cip, cjm, d, e, dip, gv);
            abs_update_s(i, j, sxx, sxy, syy, absorb_coeff);
        }
    }
}
