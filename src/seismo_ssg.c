
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

/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities or pressure or curl and div) 
 *   at receiver positions in arrays
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_ssg(int lsamp, int **recpos, float *hc, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    int i, j, ins, nxrec, nyrec, m;
    float vxx, vyy, vxy, vyx;

    float dhi = 1.0 / gv->DH;
    int fdoh = gv->FDORDER / 2;

    ins = lsamp / gv->NDT;
    for (int itr = 1; itr <= gv->NTR; itr++) {
        nxrec = recpos[1][itr];
        nyrec = recpos[2][itr];
        switch (gv->SEISMO) {
          case 1:              /* particle velocities */
              gv->SECTIONVX[itr][ins] = mpw->pvx[nyrec][nxrec];
              gv->SECTIONVY[itr][ins] = mpw->pvy[nyrec][nxrec];
              break;

          case 2:              /* pressure */
              i = nxrec;
              j = nyrec;
              //gv->SECTIONP[itr][ins]=-sxx[nyrec][nxrec]-syy[nyrec][nxrec]; // unscaled amplitude
              gv->SECTIONP[itr][ins] = ((3.0 * mpm->ppi[j][i] - 4.0 * mpm->pu[j][i]) / (2.0 * mpm->ppi[j][i] - 2.0 * mpm->pu[j][i])) * (-mpw->psxx[nyrec][nxrec] - mpw->psyy[nyrec][nxrec]) / 3;    // true amplitude
              break;

          case 3:              /* curl +div */
              i = nxrec;
              j = nyrec;

              vxx = 0;
              vyy = 0;
              vyx = 0;
              vxy = 0;
              for (m = 1; m <= fdoh; m++) {
                  vxx += hc[m] * (mpw->pvx[j][i + m - 1] - mpw->pvx[j][i - m]);
                  vyy += hc[m] * (mpw->pvy[j + m - 1][i] - mpw->pvy[j - m][i]);
                  vyx += hc[m] * (mpw->pvy[j][i + m] - mpw->pvy[j][i - m + 1]);
                  vxy += hc[m] * (mpw->pvx[j + m][i] - mpw->pvx[j - m + 1][i]);
              }
              vxx *= dhi;
              vyy *= dhi;
              vyx *= dhi;
              vxy *= dhi;

              gv->SECTIONDIV[itr][ins] = (vxx + vyy) * sqrt(mpm->ppi[j][i]);
              gv->SECTIONCURL[itr][ins] = (vxy - vyx) * sqrt(mpm->pu[j][i]);
              break;

          case 4:              /* all */
              i = nxrec;
              j = nyrec;

              vxx = 0;
              vyy = 0;
              vyx = 0;
              vxy = 0;
              for (m = 1; m <= fdoh; m++) {
                  vxx += hc[m] * (mpw->pvx[j][i + m - 1] - mpw->pvx[j][i - m]);
                  vyy += hc[m] * (mpw->pvy[j + m - 1][i] - mpw->pvy[j - m][i]);
                  vyx += hc[m] * (mpw->pvy[j][i + m] - mpw->pvy[j][i - m + 1]);
                  vxy += hc[m] * (mpw->pvx[j + m][i] - mpw->pvx[j - m + 1][i]);
              }
              vxx *= dhi;
              vyy *= dhi;
              vyx *= dhi;
              vxy *= dhi;

              gv->SECTIONDIV[itr][ins] = (vxx + vyy) * sqrt(mpm->ppi[j][i]);
              gv->SECTIONCURL[itr][ins] = (vxy - vyx) * sqrt(mpm->pu[j][i]);
              gv->SECTIONVX[itr][ins] = mpw->pvx[nyrec][nxrec];
              gv->SECTIONVY[itr][ins] = mpw->pvy[nyrec][nxrec];
              //gv->SECTIONP[itr][ins]=-sxx[nyrec][nxrec]-syy[nyrec][nxrec]; //unscaled amplitude
              gv->SECTIONP[itr][ins] = ((3.0 * mpm->ppi[j][i] - 4.0 * mpm->pu[j][i]) / (2.0 * mpm->ppi[j][i] -
                                                                                        2.0 * mpm->pu[j][i])) *
                  (-mpw->psxx[nyrec][nxrec] - mpw->psyy[nyrec][nxrec]) / 3;
              break;
        }
    }
}
