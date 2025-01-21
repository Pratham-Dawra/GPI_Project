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
 *   updating particle velocities at interior gridpoints (excluding boundarys) [GX2+1...GX3][GY2+1...GY3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_v_interior_4(int nt, float **srcpos_loc, float **signals, int nsrc, float *hc,
                         MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    float amp;
    float sxx_x, sxy_x, sxy_y, syy_y;
    double time1 = 0.0, time2 = 0.0;

    float dtdh = gv->DT / gv->DH;

    /* Coefficients for Adam Bashforth */
    float c1 = 13.0 / 12.0;
    float c2 = -5.0 / 24.0;
    float c3 = 1.0 / 6.0;
    float c4 = -1.0 / 24.0;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating particle velocities...\n");
    }

    float *svx_1_j, *svy_1_j;
    float *svx_2_j, *svy_2_j;
    float *svx_3_j, *svy_3_j;
    float *svx_4_j, *svy_4_j;

    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */

    for (int l = 1; l <= nsrc; l++) {
	int i = (int)srcpos_loc[1][l];
	int j = (int)srcpos_loc[2][l];
	float azi_rad = srcpos_loc[7][l] * PI / 180;

	//amp=signals[l][nt]; // unscaled force amplitude
	amp = (gv->DT * signals[l][nt]) / (gv->DH * gv->DH);    // scaled force amplitude with F= 1N
	
	gv->SOURCE_TYPE = (int)srcpos_loc[8][l];
	
	switch (gv->SOURCE_TYPE) {
	case 2:          /* single force in x */
	    mpw->pvx[j][i] += mpm->prip[j][i] * amp;
	    
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
	     * vx[j][i+m-1]  +=  hc[m]*mpm->prip[j][i]*amp;
	     * vx[j][i-m]    +=  hc[m]*mpm->prip[j][i-1]*amp;
	     * } */
	    break;
	case 3:          /* single force in y */
	    mpw->pvy[j][i] += mpm->prjp[j][i] * amp;
	    
	    /*for (m=1; m<=fdoh; m++) {       
	     * vy[j+m-1][i]  +=  hc[m]*mpm->prjp[j][i]*amp;
	     * vy[j-m][i]    +=  hc[m]*mpm->prjp[j][i-1]*amp;
	     * } */
	    break;
	case 4:          /* custom force */
	    mpw->pvx[j][i] += sin(azi_rad) * (mpm->prip[j][i] * amp);
	    mpw->pvy[j][i] += cos(azi_rad) * (mpm->prjp[j][i] * amp);
	    
	    /*for (m=1; m<=fdoh; m++) {
	     * vx[j][i+m-1]  +=  sin(azi_rad)*(hc[m]*mpm->prip[j][i]*amp);
	     * vx[j][i-m]    +=  sin(azi_rad)*(hc[m]*mpm->prip[j][i-1]*amp);
	     * vy[j+m-1][i]  +=  cos(azi_rad)*(hc[m]*mpm->prjp[j][i]*amp);
	     * vy[j-m][i]    +=  cos(azi_rad)*(hc[m]*mpm->prjp[j][i-1]*amp);
	     * } */
	    break;
	}
    }

    switch (gv->FDORDER) {
      case 2:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      case 4:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i])
                      + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1])
                      + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i])
                      + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i])
                      + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      case 6:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i])
                      + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1])
                      + hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1])
                      + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2])
                      + hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i])
                      + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i])
                      + hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i])
                      + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i])
                      + hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      case 8:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i])
                      + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1])
                      + hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2])
                      + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1])
                      + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2])
                      + hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3])
                      + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i])
                      + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i])
                      + hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i])
                      + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i])
                      + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i])
                      + hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i])
                      + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      case 10:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i])
                      + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1])
                      + hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2])
                      + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3])
                      + hc[5] * (mpw->psxx[j][i + 5] - mpw->psxx[j][i - 4]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1])
                      + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2])
                      + hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3])
                      + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4])
                      + hc[5] * (mpw->psxy[j][i + 4] - mpw->psxy[j][i - 5]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i])
                      + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i])
                      + hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i])
                      + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i])
                      + hc[5] * (mpw->psxy[j + 4][i] - mpw->psxy[j - 5][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i])
                      + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i])
                      + hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i])
                      + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i])
                      + hc[5] * (mpw->psyy[j + 5][i] - mpw->psyy[j - 4][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      case 12:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i])
                      + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1])
                      + hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2])
                      + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3])
                      + hc[5] * (mpw->psxx[j][i + 5] - mpw->psxx[j][i - 4])
                      + hc[6] * (mpw->psxx[j][i + 6] - mpw->psxx[j][i - 5]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1])
                      + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2])
                      + hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3])
                      + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4])
                      + hc[5] * (mpw->psxy[j][i + 4] - mpw->psxy[j][i - 5])
                      + hc[6] * (mpw->psxy[j][i + 5] - mpw->psxy[j][i - 6]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i])
                      + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i])
                      + hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i])
                      + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i])
                      + hc[5] * (mpw->psxy[j + 4][i] - mpw->psxy[j - 5][i])
                      + hc[6] * (mpw->psxy[j + 5][i] - mpw->psxy[j - 6][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i])
                      + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i])
                      + hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i])
                      + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i])
                      + hc[5] * (mpw->psyy[j + 5][i] - mpw->psyy[j - 4][i])
                      + hc[6] * (mpw->psyy[j + 6][i] - mpw->psyy[j - 5][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
      default:                 // 2nd order
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              svx_1_j = *(mpw->svx_1 + j);
              svy_1_j = *(mpw->svy_1 + j);
              svx_2_j = *(mpw->svx_2 + j);
              svy_2_j = *(mpw->svy_2 + j);
              svx_3_j = *(mpw->svx_3 + j);
              svy_3_j = *(mpw->svy_3 + j);
              svx_4_j = *(mpw->svx_4 + j);
              svy_4_j = *(mpw->svy_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]);
                  sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]);
                  sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]);
                  syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]);
                  // Save derivations
                  *(svx_1_j + i) = sxx_x + sxy_y;
                  *(svy_1_j + i) = sxy_x + syy_y;
                  mpw->pvx[j][i] +=
                      (c1 * (*(svx_1_j + i)) + c2 * (*(svx_2_j + i)) + c3 * (*(svx_3_j + i)) +
                       c4 * (*(svx_4_j + i))) * dtdh * mpm->prip[j][i];
                  mpw->pvy[j][i] +=
                      (c1 * (*(svy_1_j + i)) + c2 * (*(svy_2_j + i)) + c3 * (*(svy_3_j + i)) +
                       c4 * (*(svy_4_j + i))) * dtdh * mpm->prjp[j][i];
              }
          }
          break;
    }                           /* end of switch(gv->FDORDER) */

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating particle velocities (real time: %4.3fs).\n", time2 - time1);
    }
}
