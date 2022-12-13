
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
 *   updating particle velocities at interior gridpoints (excluding boundarys) [gx2+1...gx3][gy2+1...gy3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_v_interior(int *gx, int *gy, int nt,
                       float **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float **rip, float **rjp,
                       float **srcpos_loc, float **signals, int nsrc, float *hc, GlobVar *gv)
{
    float amp;
    float sxx_x, sxy_x, sxy_y, syy_y;
    double time1 = 0.0, time2 = 0.0;

    float dtdh = gv->DT / gv->DH;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating particle velocities...\n");
    }

    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */

    if (!gv->CHECKPTREAD)
        for (int l = 1; l <= nsrc; l++) {
            int i = (int)srcpos_loc[1][l];
            int j = (int)srcpos_loc[2][l];
            float azi_rad = srcpos_loc[7][l] * PI / 180;

            //amp=signals[l][nt]; // unscaled force amplitude
            amp = (gv->DT * signals[l][nt]) / (gv->DH * gv->DH);    // scaled force amplitude with F= 1N

            gv->SOURCE_TYPE = (int)srcpos_loc[8][l];

            switch (gv->SOURCE_TYPE) {
              case 2:          /* single force in x */

                  vx[j][i] += rip[j][i] * amp;

                  /* previous implementation of body forces as seismic sources.
                   * Implementation according to Coutant et al., BSSA, Vol. 85, No 5, 1507-1512.
                   * The stress tensor components sxx and syy are incremented prior to
                   * particle velocity update. Thereby the body force (both directions)
                   * are located at full grid point (i,j) (same position as pressure source).
                   * as a consequence, source signals are added [and weighted] at multiple grid points.
                   * This implementation works but not quite physical when considering e.g. a force of 1 N
                   * and aiming to gain the particle velocity strictly according to that force.
                   * Therefore it has been commented */

                  /*for (m=1; m<=fdoh; m++) {
                   * vx[j][i+m-1]  +=  hc[m]*rip[j][i]*amp;
                   * vx[j][i-m]    +=  hc[m]*rip[j][i-1]*amp;
                   * } */
                  break;
              case 3:          /* single force in y */
                  vy[j][i] += rjp[j][i] * amp;
                  /*for (m=1; m<=fdoh; m++) {
                   * vy[j+m-1][i]  +=  hc[m]*rjp[j][i]*amp;
                   * vy[j-m][i]    +=  hc[m]*rjp[j][i-1]*amp;
                   * } */
                  break;
              case 4:          /* custom force */
                  vx[j][i] += sin(azi_rad) * (rip[j][i] * amp);
                  vy[j][i] += cos(azi_rad) * (rjp[j][i] * amp);
                  /*for (m=1; m<=fdoh; m++) {
                   * vx[j][i+m-1]  +=  sin(azi_rad)*(hc[m]*rip[j][i]*amp);
                   * vx[j][i-m]    +=  sin(azi_rad)*(hc[m]*rip[j][i-1]*amp);
                   * vy[j+m-1][i]  +=  cos(azi_rad)*(hc[m]*rjp[j][i]*amp);
                   * vy[j-m][i]    +=  cos(azi_rad)*(hc[m]*rjp[j][i-1]*amp);
                   * } */
                  break;
            }
        }

    switch (gv->FDORDER) {
      case 2:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      case 4:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      case 6:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i])
                      + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1])
                      + hc[3] * (sxx[j][i + 3] - sxx[j][i - 2]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1])
                      + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2])
                      + hc[3] * (sxy[j][i + 2] - sxy[j][i - 3]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i])
                      + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i])
                      + hc[3] * (sxy[j + 2][i] - sxy[j - 3][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i])
                      + hc[2] * (syy[j + 2][i] - syy[j - 1][i])
                      + hc[3] * (syy[j + 3][i] - syy[j - 2][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      case 8:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i])
                      + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1])
                      + hc[3] * (sxx[j][i + 3] - sxx[j][i - 2])
                      + hc[4] * (sxx[j][i + 4] - sxx[j][i - 3]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1])
                      + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2])
                      + hc[3] * (sxy[j][i + 2] - sxy[j][i - 3])
                      + hc[4] * (sxy[j][i + 3] - sxy[j][i - 4]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i])
                      + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i])
                      + hc[3] * (sxy[j + 2][i] - sxy[j - 3][i])
                      + hc[4] * (sxy[j + 3][i] - sxy[j - 4][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i])
                      + hc[2] * (syy[j + 2][i] - syy[j - 1][i])
                      + hc[3] * (syy[j + 3][i] - syy[j - 2][i])
                      + hc[4] * (syy[j + 4][i] - syy[j - 3][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      case 10:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i])
                      + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1])
                      + hc[3] * (sxx[j][i + 3] - sxx[j][i - 2])
                      + hc[4] * (sxx[j][i + 4] - sxx[j][i - 3])
                      + hc[5] * (sxx[j][i + 5] - sxx[j][i - 4]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1])
                      + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2])
                      + hc[3] * (sxy[j][i + 2] - sxy[j][i - 3])
                      + hc[4] * (sxy[j][i + 3] - sxy[j][i - 4])
                      + hc[5] * (sxy[j][i + 4] - sxy[j][i - 5]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i])
                      + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i])
                      + hc[3] * (sxy[j + 2][i] - sxy[j - 3][i])
                      + hc[4] * (sxy[j + 3][i] - sxy[j - 4][i])
                      + hc[5] * (sxy[j + 4][i] - sxy[j - 5][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i])
                      + hc[2] * (syy[j + 2][i] - syy[j - 1][i])
                      + hc[3] * (syy[j + 3][i] - syy[j - 2][i])
                      + hc[4] * (syy[j + 4][i] - syy[j - 3][i])
                      + hc[5] * (syy[j + 5][i] - syy[j - 4][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      case 12:
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i])
                      + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1])
                      + hc[3] * (sxx[j][i + 3] - sxx[j][i - 2])
                      + hc[4] * (sxx[j][i + 4] - sxx[j][i - 3])
                      + hc[5] * (sxx[j][i + 5] - sxx[j][i - 4])
                      + hc[6] * (sxx[j][i + 6] - sxx[j][i - 5]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1])
                      + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2])
                      + hc[3] * (sxy[j][i + 2] - sxy[j][i - 3])
                      + hc[4] * (sxy[j][i + 3] - sxy[j][i - 4])
                      + hc[5] * (sxy[j][i + 4] - sxy[j][i - 5])
                      + hc[6] * (sxy[j][i + 5] - sxy[j][i - 6]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i])
                      + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i])
                      + hc[3] * (sxy[j + 2][i] - sxy[j - 3][i])
                      + hc[4] * (sxy[j + 3][i] - sxy[j - 4][i])
                      + hc[5] * (sxy[j + 4][i] - sxy[j - 5][i])
                      + hc[6] * (sxy[j + 5][i] - sxy[j - 6][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i])
                      + hc[2] * (syy[j + 2][i] - syy[j - 1][i])
                      + hc[3] * (syy[j + 3][i] - syy[j - 2][i])
                      + hc[4] * (syy[j + 4][i] - syy[j - 3][i])
                      + hc[5] * (syy[j + 5][i] - syy[j - 4][i])
                      + hc[6] * (syy[j + 6][i] - syy[j - 5][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
      default:                 // 2nd order
          for (int j = gy[2] + 1; j <= gy[3]; j++) {
              for (int i = gx[2] + 1; i <= gx[3]; i++) {
                  sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i]);
                  sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1]);
                  sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i]);
                  syy_y = hc[1] * (syy[j + 1][i] - syy[j][i]);
                  vx[j][i] += (sxx_x + sxy_y) * dtdh * rip[j][i];
                  vy[j][i] += (sxy_x + syy_y) * dtdh * rjp[j][i];
              }
          }
          break;
    }                           /* end of switch(gv->FDORDER) */

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating particle velocities (real time: %4.3fs).\n", time2 - time1);
    }
}
