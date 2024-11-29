
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
 * store amplitudes (particle velocities or pressure or curl and div)
 * at receiver positions in arrays
 *----------------------------------------------------------------------*/

#include "fd.h"

void seismo_ssg(int lsamp, int **recpos, MemModel *mpm, MemWavefield *mpw, GlobVar *gv)
{
    int ins, nxrec, nyrec;

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
              gv->SECTIONP[itr][ins] = -mpw->psxx[nyrec][nxrec] - mpw->psyy[nyrec][nxrec]; // unscaled amplitude
              /*gv->SECTIONP[itr][ins] = ((3.0 * mpm->ppi[j][i] - 4.0 * mpm->pu[j][i]) / (2.0 * mpm->ppi[j][i] - 2.0 * mpm->pu[j][i]))
                                       * (-mpw->psxx[nyrec][nxrec] - mpw->psyy[nyrec][nxrec]) / 3;  */  // true amplitude
              break;

          case 3:              /* curl +div */

              gv->SECTIONDIV[itr][ins] = (mpw->pvxx[nyrec][nxrec] + mpw->pvyy[nyrec][nxrec]) * sqrt(mpm->ppi[nyrec][nxrec]);
              gv->SECTIONCURL[itr][ins] = (mpw->pvxy[nyrec][nxrec] - mpw->pvyx[nyrec][nxrec]) * sqrt(mpm->pu[nyrec][nxrec]);
              break;

          case 4:              /* all */

              gv->SECTIONDIV[itr][ins] = (mpw->pvxx[nyrec][nxrec] + mpw->pvyy[nyrec][nxrec]) * sqrt(mpm->ppi[nyrec][nxrec]);
              gv->SECTIONCURL[itr][ins] = (mpw->pvxy[nyrec][nxrec] - mpw->pvyx[nyrec][nxrec]) * sqrt(mpm->pu[nyrec][nxrec]);
              gv->SECTIONVX[itr][ins] = mpw->pvx[nyrec][nxrec];
              gv->SECTIONVY[itr][ins] = mpw->pvy[nyrec][nxrec];
              //gv->SECTIONP[itr][ins]=-sxx[nyrec][nxrec]-syy[nyrec][nxrec]; //unscaled amplitude
              gv->SECTIONP[itr][ins] = ((3.0 * mpm->ppi[nyrec][nxrec] - 4.0 * mpm->pu[nyrec][nxrec]) / (2.0 * mpm->ppi[nyrec][nxrec] -
                                         2.0 * mpm->pu[nyrec][nxrec])) * (-mpw->psxx[nyrec][nxrec] - mpw->psyy[nyrec][nxrec]) / 3;
              break;
        }
    }
}
