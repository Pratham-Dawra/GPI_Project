
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
 *   updating stress components at interior gridpoints (excluding boundarys) [gv->GX2+1...gv->GX3][gv->GY2+1...gv->GY3]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_s_elastic_interior(int nt, float *hc, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    float fipjp, f, g;
    float vxx, vyy, vxy, vyx;
    double time1 = 0.0, time2 = 0.0;

    float dhi = 1.0f / gv->DH;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating stress components...\n");
    }

    switch (gv->FDORDER) {
      case 2:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  vxx = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) * dhi;
                  vyy = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) * dhi;
                  vyx = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) * dhi;
                  vxy = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) * dhi;

                  /* updating components of the stress tensor, partially */
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      case 4:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  vxx = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                         + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])) * dhi;
                  vyy = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                         + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])) * dhi;
                  vyx = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                         + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])) * dhi;
                  vxy = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                         + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])) * dhi;
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      case 6:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  vxx = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                         + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                         + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])) * dhi;
                  vyy = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                         + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                         + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])) * dhi;
                  vyx = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                         + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                         + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])) * dhi;
                  vxy = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                         + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                         + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])) * dhi;
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      case 8:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  vxx = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                         + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                         + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                         + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])) * dhi;
                  vyy = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                         + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                         + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                         + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])) * dhi;
                  vyx = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                         + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                         + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                         + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])) * dhi;
                  vxy = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                         + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                         + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                         + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])) * dhi;
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      case 10:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  vxx = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                         + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                         + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                         + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])
                         + hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5])) * dhi;
                  vyy = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                         + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                         + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                         + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])
                         + hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i])) * dhi;
                  vyx = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                         + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                         + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                         + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])
                         + hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4])) * dhi;
                  vxy = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                         + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                         + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                         + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])
                         + hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i])) * dhi;
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      case 12:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  vxx = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                         + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                         + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                         + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])
                         + hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5])
                         + hc[6] * (mpw->pvx[j][i + 5] - mpw->pvx[j][i - 6])) * dhi;
                  vyy = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                         + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                         + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                         + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])
                         + hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i])
                         + hc[6] * (mpw->pvy[j + 5][i] - mpw->pvy[j - 6][i])) * dhi;
                  vyx = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                         + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                         + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                         + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])
                         + hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4])
                         + hc[6] * (mpw->pvy[j][i + 6] - mpw->pvy[j][i - 5])) * dhi;
                  vxy = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                         + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                         + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                         + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])
                         + hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i])
                         + hc[6] * (mpw->pvx[j + 6][i] - mpw->pvx[j - 5][i])) * dhi;

                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;
                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
      default:                 /* CASE 2 */
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {

                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  vxx = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) * dhi;
                  vyy = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) * dhi;
                  vyx = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) * dhi;
                  vxy = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) * dhi;

                  /* updating components of the stress tensor, partially */
                  fipjp = mpm->puipjp[j][i] * gv->DT;
                  f = mpm->pu[j][i] * gv->DT;
                  g = mpm->ppi[j][i] * gv->DT;

                  mpw->psxy[j][i] += (fipjp * (vxy + vyx));
                  mpw->psxx[j][i] += (g * (vxx + vyy)) - (2.0 * f * vyy);
                  mpw->psyy[j][i] += (g * (vxx + vyy)) - (2.0 * f * vxx);
              }
          }
          break;
    }                           /* end of switch(gv->FDORDER) */

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating stress components (real time: %4.3fs).\n", time2 - time1);
    }
}
