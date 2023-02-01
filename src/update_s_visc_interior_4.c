
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
 *   updating stress components at interior gridpoints (excluding boundarys) [GX2+1...GX3][GY2+1...GY3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_s_visc_interior_4(int nt, float *hc, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    float sumr = 0.0f, sump = 0.0f, sumq = 0.0f;
    float ctemp;
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
    float **r_j, **p_j, **q_j, **r_2_j, **p_2_j, **q_2_j, **r_3_j, **p_3_j, **q_3_j, **r_4_j, **p_4_j, **q_4_j;
    float *r_ji, *p_ji, *q_ji, *r_2_ji, *p_2_ji, *q_2_ji, *r_3_ji, *p_3_ji, *q_3_ji, *r_4_ji, *p_4_ji, *q_4_ji;

    float sumxx = 0.0f, sumyy = 0.0f, sumxy = 0.0f, sumyx = 0.0f;

    float dthalbe = gv->DT / 2.0;
    float dhi = 1.0 / gv->DH;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating stress components...\n");
    }

    switch (gv->FDORDER) {      /* standard staggered grid (SSG) */
      case 2:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
                  /* Compute values for shearmodulus u[j][i],
                   * P-wave modulus pi[j][i],
                   * tau for S-waves and P-waves taus[j][i],
                   * taup[j][i] at staggered grid points: */

                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  *(vxx_1_j + i) = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]);
                  *(vyy_1_j + i) = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]);
                  *(vyx_1_j + i) = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]);
                  *(vxy_1_j + i) = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] = mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) +
                                                                                                      (*(r_2_ji + l))) -
                                                          ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                                          ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                                          (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;

      case 4:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);

                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]));
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]));
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]));
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]));

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;

      case 6:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
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

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;

      case 8:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
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

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;
      case 10:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
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

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;
      case 12:
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
                  *(vxx_1_j + i) = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1])
                                    + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])
                                    + hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])
                                    + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])
                                    + hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5])
                                    + hc[6] * (mpw->pvx[j][i + 5] - mpw->pvx[j][i - 6])
                      );
                  *(vyy_1_j + i) = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i])
                                    + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])
                                    + hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])
                                    + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])
                                    + hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i])
                                    + hc[6] * (mpw->pvy[j + 5][i] - mpw->pvy[j - 6][i])
                      );
                  *(vyx_1_j + i) = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i])
                                    + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])
                                    + hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])
                                    + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])
                                    + hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4])
                                    + hc[6] * (mpw->pvy[j][i + 6] - mpw->pvy[j][i - 5])
                      );
                  *(vxy_1_j + i) = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i])
                                    + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])
                                    + hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])
                                    + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])
                                    + hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i])
                                    + hc[6] * (mpw->pvx[j + 6][i] - mpw->pvx[j - 5][i])
                      );

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;
      default:                 /* Case 2 */
          for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
              r_j = *(mpw->pr + j);
              r_2_j = *(mpw->pr_2 + j);
              r_3_j = *(mpw->pr_3 + j);
              r_4_j = *(mpw->pr_4 + j);
              q_j = *(mpw->pq + j);
              q_2_j = *(mpw->pq_2 + j);
              q_3_j = *(mpw->pq_3 + j);
              q_4_j = *(mpw->pq_4 + j);
              p_j = *(mpw->pp + j);
              p_2_j = *(mpw->pp_2 + j);
              p_3_j = *(mpw->pp_3 + j);
              p_4_j = *(mpw->pp_4 + j);
              vxx_1_j = *(mpw->vxx_1 + j), vyy_1_j = *(mpw->vyy_1 + j), vxy_1_j = *(mpw->vxy_1 + j), vyx_1_j =
                  *(mpw->vyx_1 + j);
              vxx_2_j = *(mpw->vxx_2 + j), vyy_2_j = *(mpw->vyy_2 + j), vxy_2_j = *(mpw->vxy_2 + j), vyx_2_j =
                  *(mpw->vyx_2 + j);
              vxx_3_j = *(mpw->vxx_3 + j), vyy_3_j = *(mpw->vyy_3 + j), vxy_3_j = *(mpw->vxy_3 + j), vyx_3_j =
                  *(mpw->vyx_3 + j);
              vxx_4_j = *(mpw->vxx_4 + j), vyy_4_j = *(mpw->vyy_4 + j), vxy_4_j = *(mpw->vxy_4 + j), vyx_4_j =
                  *(mpw->vyx_4 + j);
              for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                  r_ji = *(r_j + i);
                  r_2_ji = *(r_2_j + i);
                  r_3_ji = *(r_3_j + i);
                  r_4_ji = *(r_4_j + i);
                  q_ji = *(q_j + i);
                  q_2_ji = *(q_2_j + i);
                  q_3_ji = *(q_3_j + i);
                  q_4_ji = *(q_4_j + i);
                  p_ji = *(p_j + i);
                  p_2_ji = *(p_2_j + i);
                  p_3_ji = *(p_3_j + i);
                  p_4_ji = *(p_4_j + i);
                  /* Compute values for shearmodulus u[j][i],
                   * P-wave modulus pi[j][i],
                   * tau for S-waves and P-waves taus[j][i],
                   * taup[j][i] at staggered grid points: */

                  /* spatial derivatives of the components of the velocities */
                  /* using Holberg coefficients */
                  *(vxx_1_j + i) = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]);
                  *(vyy_1_j + i) = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]);
                  *(vyx_1_j + i) = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]);
                  *(vxy_1_j + i) = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      sumr += c1 * (*(r_ji + l)) + c2 * (*(r_2_ji + l)) + c3 * (*(r_3_ji + l)) + c4 * (*(r_4_ji + l));
                      sump += c1 * (*(p_ji + l)) + c2 * (*(p_2_ji + l)) + c3 * (*(p_3_ji + l)) + c4 * (*(p_4_ji + l));
                      sumq += c1 * (*(q_ji + l)) + c2 * (*(q_2_ji + l)) + c3 * (*(q_3_ji + l)) + c4 * (*(q_4_ji + l));
                  }

                  // Calculate Adams-Bashforth stuff
                  sumxx = c1 * (*(vxx_1_j + i)) + c2 * (*(vxx_2_j + i)) + c3 * (*(vxx_3_j + i)) + c4 * (*(vxx_4_j + i));
                  sumyy = c1 * (*(vyy_1_j + i)) + c2 * (*(vyy_2_j + i)) + c3 * (*(vyy_3_j + i)) + c4 * (*(vyy_4_j + i));
                  sumxy = c1 * (*(vxy_1_j + i)) + c2 * (*(vxy_2_j + i)) + c3 * (*(vxy_3_j + i)) + c4 * (*(vxy_4_j + i));
                  sumyx = c1 * (*(vyx_1_j + i)) + c2 * (*(vyx_2_j + i)) + c3 * (*(vyx_3_j + i)) + c4 * (*(vyx_4_j + i));

                  /* updating components of the stress tensor, partially */
                  mpw->psxy[j][i] += (mpm->fipjp[j][i] * (sumxy + sumyx)) * dhi + (dthalbe * sumr);
                  mpw->psxx[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumyy * dhi) + (dthalbe * sump);
                  mpw->psyy[j][i] +=
                      (mpm->g[j][i] * (sumxx + sumyy) * dhi) - (2.0 * mpm->f[j][i] * sumxx * dhi) + (dthalbe * sumq);

                  sumr = sump = sumq = 0.0f;
                  for (int l = 1; l <= gv->L; l++) {
                      ctemp = 2 * (1 - mpm->cip[l]) / c1;
                      mpw->pr_4[j][i][l] =
                          mpm->bip[l] * ((*(r_ji + l)) * mpm->cip[l] - ctemp * c2 * ((*(r_ji + l)) + (*(r_2_ji + l))) -
                                         ctemp * c3 * ((*(r_2_ji + l)) + (*(r_3_ji + l))) -
                                         ctemp * c4 * ((*(r_3_ji + l)) + (*(r_4_ji + l))) -
                                         (mpm->dip[j][i][l] * (sumxy + sumyx) * dhi));
                      mpw->pp_4[j][i][l] =
                          mpm->bjm[l] * ((*(p_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(p_ji + l)) + (*(p_2_ji + l))) -
                                         ctemp * c3 * ((*(p_2_ji + l)) + (*(p_3_ji + l))) -
                                         ctemp * c4 * ((*(p_3_ji + l)) + (*(p_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumyy * dhi));
                      mpw->pq_4[j][i][l] =
                          mpm->bjm[l] * ((*(q_ji + l)) * mpm->cjm[l] - ctemp * c2 * ((*(q_ji + l)) + (*(q_2_ji + l))) -
                                         ctemp * c3 * ((*(q_2_ji + l)) + (*(q_3_ji + l))) -
                                         ctemp * c4 * ((*(q_3_ji + l)) + (*(q_4_ji + l))) -
                                         (mpm->e[j][i][l] * (sumxx + sumyy) * dhi) +
                                         (2.0 * mpm->d[j][i][l] * sumxx * dhi));
                      sumr += c1 * (*(r_4_ji + l)) + c2 * (*(r_ji + l)) + c3 * (*(r_2_ji + l)) + c4 * (*(r_3_ji + l));
                      sump += c1 * (*(p_4_ji + l)) + c2 * (*(p_ji + l)) + c3 * (*(p_2_ji + l)) + c4 * (*(p_3_ji + l));
                      sumq += c1 * (*(q_4_ji + l)) + c2 * (*(q_ji + l)) + c3 * (*(q_2_ji + l)) + c4 * (*(q_3_ji + l));
                  }

                  mpw->psxy[j][i] += (dthalbe * sumr);
                  mpw->psxx[j][i] += (dthalbe * sump);
                  mpw->psyy[j][i] += (dthalbe * sumq);
              }
          }
          break;
    }                           /* end of switch(gv->FDORDER) */

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating stress components (real time: %4.3fs).\n", time2 - time1);
    }
}
