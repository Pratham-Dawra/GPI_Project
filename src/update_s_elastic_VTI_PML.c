
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
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_VTI_PML(MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    int h1;
    float vxx, vyy, vxy, vyx;

    /* left boundary */
    for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            cpml_update_s_x(i, j, &vxx, &vyx, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* right boundary */
    for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            h1 = (i - gv->NX + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* top boundary */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            cpml_update_s_y(i, j, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* bottom boundary */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            h1 = (j - gv->NY + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* corners */

    /*left-top */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            cpml_update_s_x(i, j, &vxx, &vyx, mpm, mpw);
            cpml_update_s_y(i, j, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /*left-bottom */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            cpml_update_s_x(i, j, &vxx, &vyx, mpm, mpw);
            h1 = (j - gv->NY + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* right-top */
    for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            h1 = (i - gv->NX + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, mpm, mpw);
            cpml_update_s_y(i, j, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }

    /* right-bottom */
    for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
        for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
            gv->FDOP_S(i, j, &vxx, &vyx, &vxy, &vyy, mpw);
            h1 = (i - gv->NX + 2 * gv->FW);
            cpml_update_s_x(h1, j, &vxx, &vyx, mpm, mpw);
            h1 = (j - gv->NY + 2 * gv->FW);
            cpml_update_s_y(i, h1, &vxy, &vyy, mpm, mpw);
            wavefield_update_s_el_vti(i, j, vxx, vyx, vxy, vyy, mpm, mpw);
        }
    }
}
