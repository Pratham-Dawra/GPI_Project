
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

/*----------------------------------------------------------------------
 * Write 2D snapshot for current timestep to file
 *----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void snap(int nt, int nsnap, float *hc, MemModel *mpm, MemWavefield *mpw, GlobVar *gv)
{
    /* 
     * different data formats of output:
     * SNAP_FORMAT=1  :  SU (IEEE)
     * SNAP_FORMAT=2  :  ASCII
     * SNAP_FORMAT=3  :  BINARY (IEEE)
     * 
     * different types:
     * SNAP=1 : values in vx and vy
     * SNAP=2 : -(vx+vy) (pressure field)
     * SNAP=3 : divergence of vx and vy (energy of compressional waves)
     * and curl of vx and vy (energy of shear waves)
     * SNAP=4 : both particle velocities (type=1) and energy (type=3)
     */

    float amp, vyx, vxy, vxx, vyy;
    float **divfield, **curlfield;
    char snapfile_x[STRING_SIZE * 2], snapfile_y[STRING_SIZE * 2], snapfile_div[STRING_SIZE * 2];
    char snapfile_rot[STRING_SIZE * 2], snapfile_p[STRING_SIZE * 2], ext[8];
    FILE *fpx1, *fpy1, *fpx2, *fpy2, *fpp;

    float dhi = 1.0 / gv->DH;
    int fdoh = gv->FDORDER / 2;

    switch (gv->SNAP_FORMAT) {
      case 1:
          sprintf(ext, ".su");
          break;
      case 2:
          sprintf(ext, ".asc");
          break;
      case 3:
          sprintf(ext, ".bin");
          break;
    }
    sprintf(snapfile_x, "%s%s.vx.%i.%i", gv->SNAP_FILE, ext, gv->POS[1], gv->POS[2]);
    sprintf(snapfile_y, "%s%s.vy.%i.%i", gv->SNAP_FILE, ext, gv->POS[1], gv->POS[2]);
    sprintf(snapfile_div, "%s%s.div.%i.%i", gv->SNAP_FILE, ext, gv->POS[1], gv->POS[2]);
    sprintf(snapfile_rot, "%s%s.curl.%i.%i", gv->SNAP_FILE, ext, gv->POS[1], gv->POS[2]);
    sprintf(snapfile_p, "%s%s.p.%i.%i", gv->SNAP_FILE, ext, gv->POS[1], gv->POS[2]);

    log_debug("Writing snapshot data at T=%fs, file suffix '%d.%d'.\n", nt * gv->DT, gv->POS[1], gv->POS[2]);

    /* do NOT remove the FALLTHRU comment below, it is used to tell the compiler
     * that this is an intentional fall through */
    switch (gv->SNAP) {
      case 1:
          //log_debug("%s\n", snapfile_x);
          //log_debug("%s\n", snapfile_y);      
          if (nsnap == 1) {
              fpx1 = fopen(snapfile_x, "w");
              fpy1 = fopen(snapfile_y, "w");
          } else {
              fpx1 = fopen(snapfile_x, "a");
              fpy1 = fopen(snapfile_y, "a");
          }
          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {
                  writedsk(fpx1, mpw->pvx[j][i], gv->SNAP_FORMAT);
                  writedsk(fpy1, mpw->pvy[j][i], gv->SNAP_FORMAT);
              }
          }
          fclose(fpx1);
          fclose(fpy1);
          break;
      case 2:
          //log_debug("%s\n",snapfile_p);
          if (nsnap == 1) {
              fpx1 = fopen(snapfile_p, "w");
          } else {
              fpx1 = fopen(snapfile_p, "a");
          }

          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {
                  amp = -mpw->psxx[j][i] - mpw->psyy[j][i];
                  writedsk(fpx1, amp, gv->SNAP_FORMAT);
              }
          }
          fclose(fpx1);
          break;
      case 4:
          //log_debug("%s\n", snapfile_x);
          //log_debug("%s\n", snapfile_y);
          //log_debug("%s\n",snapfile_p);
          if (nsnap == 1) {
              fpx1 = fopen(snapfile_x, "w");
              fpy1 = fopen(snapfile_y, "w");
              fpp = fopen(snapfile_p, "w");
          } else {
              fpx1 = fopen(snapfile_x, "a");
              fpy1 = fopen(snapfile_y, "a");
              fpp = fopen(snapfile_p, "a");
          }

          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {
                  writedsk(fpx1, mpw->pvx[j][i], gv->SNAP_FORMAT);
                  writedsk(fpy1, mpw->pvy[j][i], gv->SNAP_FORMAT);
                  amp = -mpw->psxx[j][i] - mpw->psyy[j][i];
                  writedsk(fpp, amp, gv->SNAP_FORMAT);
              }
          }
          fclose(fpx1);
          fclose(fpy1);
          fclose(fpp);
          /* FALLTHRU */
      case 3:
          /* output of the curl of the velocity field according to Dougherty and
           * Stephen (PAGEOPH, 1988) */
          /*if (gv->NY1<=2) error("NY1 must be greater than 2."); */
          //log_debug("%s\n", snapfile_div);
          //log_debug("%s\n", snapfile_rot);
          if (nsnap == 1) {
              fpx2 = fopen(snapfile_div, "w");
              fpy2 = fopen(snapfile_rot, "w");
          } else {
              fpx2 = fopen(snapfile_div, "a");
              fpy2 = fopen(snapfile_rot, "a");
          }

          int nd = gv->FDORDER / 2;
          curlfield = matrix(-nd + 1, gv->NY + nd, -nd + 1, gv->NX + nd);

          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {

                  /* spatial derivatives using Holberg coefficients */
                  vyx = 0;
                  vxy = 0;
                  for (int m = 1; m <= fdoh; m++) {
                      vyx += hc[m] * (mpw->pvy[j][i + m] - mpw->pvy[j][i - m + 1]);
                      vxy += hc[m] * (mpw->pvx[j + m][i] - mpw->pvx[j - m + 1][i]);
                  }
                  vyx *= dhi;
                  vxy *= dhi;

                  curlfield[j][i] = (vxy - vyx) * sqrt(mpm->pu[j][i]);
              }
          }

          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {
                  writedsk(fpy2, curlfield[j][i], gv->SNAP_FORMAT);
              }
          }
          free_matrix(curlfield, -nd + 1, gv->NY + nd, -nd + 1, gv->NX + nd);

          /* output of the divergence of the velocity field according to Dougherty and
           * Stephen (PAGEOPH, 1988) */
          divfield = matrix(-nd + 1, gv->NY + nd, -nd + 1, gv->NX + nd);
          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {

                  /* spatial derivatives using Holberg coefficients */
                  vxx = 0;
                  vyy = 0;
                  for (int m = 1; m <= fdoh; m++) {
                      vxx += hc[m] * (mpw->pvx[j][i + m - 1] - mpw->pvx[j][i - m]);
                      vyy += hc[m] * (mpw->pvy[j + m - 1][i] - mpw->pvy[j - m][i]);
                  }
                  vxx *= dhi;
                  vyy *= dhi;

                  divfield[j][i] = (vxx + vyy) * sqrt(mpm->ppi[j][i]);
              }
          }

          for (int i = gv->SNAPIDX[1]; i <= gv->SNAPIDX[2]; i += gv->IDX) {
              for (int j = gv->SNAPIDX[3]; j <= gv->SNAPIDX[4]; j += gv->IDY) {
                  writedsk(fpx2, divfield[j][i], gv->SNAP_FORMAT);
              }
          }

          free_matrix(divfield, -nd + 1, gv->NY + nd, -nd + 1, gv->NX + nd);
          fclose(fpx2);
          fclose(fpy2);
          break;
    }

}
