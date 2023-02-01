
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
 * along with SOFI2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating stress components at interior gridpoints (excluding boundarys) [GX2+1...GX3][GY2+1...GY3]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_s_elastic_interior_4(int nt, float *hc, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    float fipjp, f, g;
    double time1 = 0.0, time2 = 0.0;
    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;

    float *vxx_1_j, *vyy_1_j, *vxy_1_j, *vyx_1_j;
    float *vxx_2_j, *vyy_2_j, *vxy_2_j, *vyx_2_j;
    float *vxx_3_j, *vyy_3_j, *vxy_3_j, *vyx_3_j;
    float *vxx_4_j, *vyy_4_j, *vxy_4_j, *vyx_4_j;

    float sumxx = 0.0, sumyy = 0.0, sumxy = 0.0, sumyx = 0.0;

    float dhi = 1.0 / gv->DH;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating stress components...\n");
    }

    switch (gv->FDORDER) {
      case 2:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  *(vxx_1_j + i) = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]);
                  *(vyy_1_j + i) = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]);
                  *(vyx_1_j + i) = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]);
                  *(vxy_1_j + i) = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]);
                  /* updating components of the stress tensor, partially */
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));
                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      case 4:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]));
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      case 6:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                                    + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                                    + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                                    + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                                    + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i]));
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      case 8:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                                    + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                                    + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                                    + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                                    + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                                    + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                                    + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                                    + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                                    + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i]));
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      case 10:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                                    + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                                    + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])
                                    + hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                                    + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                                    + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])
                                    + hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                                    + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                                    + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])
                                    + hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                                    + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                                    + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])
                                    + hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i]));
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      case 12:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                                    + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                                    + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])
                                    + hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5])
                                    + hc[6] * (mpw->pvx[j][i + 5] - mpw->pvx[j][i - 6]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                                    + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                                    + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])
                                    + hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i])
                                    + hc[6] * (mpw->pvy[j + 5][i] - mpw->pvy[j - 6][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                                    + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                                    + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])
                                    + hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4])
                                    + hc[6] * (mpw->pvy[j][i + 6] - mpw->pvy[j][i - 5]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                                    + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                                    + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])
                                    + hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i])
                                    + hc[6] * (mpw->pvx[j + 6][i] - mpw->pvx[j - 5][i]));
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
      default:                 /* CASE 2 */
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  *(vxx_1_j + i) = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]);
                  *(vyy_1_j + i) = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]);
                  *(vyx_1_j + i) = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]);
                  *(vxy_1_j + i) = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]);

                  /* updating components of the stress tensor, partially */
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  // Update stress
                  mpw->psxy[j][i] += (fipjp * (sumxy + sumyx)) * dhi;
                  mpw->psxx[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumyy)) * dhi;
                  mpw->psyy[j][i] += ((g * (sumxx + sumyy)) - (2.0 * f * sumxx)) * dhi;
              }
          }
          break;
    }                           /* end of switch(gv->FDORDER) */

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating stress components (real time: %4.3fs).\n", time2 - time1);
    }
}
