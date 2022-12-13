
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_VTI_PML(int nx2, int ny2, int *gx, int *gy,
                           float **vx, float **vy, float **sxx, float **syy, float **sxy,
                           float ***pr, float ***pp, float ***pq,
                           float **pc55ipjpu, float **pc13u, float **pc11u, float **pc33u,
                           float ***pc55ipjpd, float ***pc13d, float ***pc11d, float ***pc33d,
                           float *bip, float *cip,
                           float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                           float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                           float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx, GlobVar *gv)
{
    int h1;
    float vxx, vyy, vxy, vyx;

    /* left boundary */
    for (int j = gy[2] + 1; j <= gy[3]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            cpml_update_s_x(i, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* right boundary */
    for (int j = gy[2] + 1; j <= gy[3]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* top boundary */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[2] + 1; i <= gx[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            cpml_update_s_y(i, j, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* bottom boundary */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[2] + 1; i <= gx[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* corners */

    /*left-top */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            cpml_update_s_x(i, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            cpml_update_s_y(i, j, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /*left-bottom */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[1]; i <= gx[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            cpml_update_s_x(i, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* right-top */
    for (int j = gy[1]; j <= gy[2]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            cpml_update_s_y(i, j, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }

    /* right-bottom */
    for (int j = gy[3] + 1; j <= gy[4]; j++) {
        for (int i = gx[3] + 1; i <= gx[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, vx, vy);
            h1 = (i - nx2 + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, psi_vxx, psi_vyx);
            h1 = (j - ny2 + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vxy);
            wavefield_update_s_visc_VTI(i, j, vxx, vyx, vxy, vyy, sxy, sxx, syy, pr, pp, pq,
                                        pc55ipjpu, pc13u, pc11u, pc33u, pc55ipjpd, pc13d, pc11d, pc33d, bip, cip, gv);
        }
    }
}
