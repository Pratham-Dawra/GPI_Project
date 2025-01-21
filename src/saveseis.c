
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

/* ------------------------------------------------------------------------
 *   Write seismograms to file(s)
 * ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis(int ishot, AcqVar *acq, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{

    /* write seismograms to file(s) */
    if (gv->SEISMO) {

        /* saves seismograms portion of each PE individually to file */
        //if (gv->NTR> 0) saveseis(SECTIONVX,SECTIONVY,SECTIONP,SECTIONCURL,SECTIONDIV,acq->recpos,acq->recpos_loc,gv->NTR,acq->srcpos_current,ishot,gv->NS);

        /* merge of seismogram data from all PE and output data collectively */
        switch (gv->SEISMO) {
          case 1:              /* particle velocities only */
              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionvxcalc[i][j] = gv->SECTIONVX[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONVX, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  if (gv->MODE == FWI && vinv->LNORM == 8) {    /* Why only for SEISMO==1 ??? */
                      calc_envelope(gv->SEISMO_FULLDATA, gv->SEISMO_FULLDATA, gv->NTRG, gv->NS);
                  }
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 1, gv);
              }

              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionvycalc[i][j] = gv->SECTIONVY[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONVY, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  if (gv->MODE == FWI && vinv->LNORM == 8) {    /* Why only for SEISMO==1 ??? */
                      calc_envelope(gv->SEISMO_FULLDATA, gv->SEISMO_FULLDATA, gv->NTRG, gv->NS);
                  }
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 2, gv);
              }

              break;
          case 2:              /* pressure only */
              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionpcalc[i][j] = gv->SECTIONP[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONP, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 4, gv);
              }

              break;
          case 3:              /* curl and div only */
              catseis(gv->SECTIONDIV, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0)
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 5, gv);
              catseis(gv->SECTIONCURL, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0)
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 6, gv);

              break;
          case 4:              /* everything */
              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionvxcalc[i][j] = gv->SECTIONVX[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONVX, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  if (gv->MODE == FWI && vinv->LNORM == 8) {    /* Why only for SEISMO==1 ??? */
                      calc_envelope(gv->SEISMO_FULLDATA, gv->SEISMO_FULLDATA, gv->NTRG, gv->NS);
                  }
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 1, gv);
              }

              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionvycalc[i][j] = gv->SECTIONVY[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONVY, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  if (gv->MODE == FWI && vinv->LNORM == 8) {    /* Why only for SEISMO==1 ??? */
                      calc_envelope(gv->SEISMO_FULLDATA, gv->SEISMO_FULLDATA, gv->NTRG, gv->NS);
                  }
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 2, gv);
              }

              if (gv->MODE == FWI) {
                  for (int i = 1; i <= gv->NTR; i++) {
                      for (int j = 1; j <= gv->NS; j++) {
                          minv->sectionpcalc[i][j] = gv->SECTIONP[i][j];
                      }
                  }
              }
              catseis(gv->SECTIONP, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0) {
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 4, gv);
              }

              catseis(gv->SECTIONDIV, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0)
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 5, gv);

              catseis(gv->SECTIONCURL, gv->SEISMO_FULLDATA, acq->recswitch, gv->NTRG, gv->NS);
              if (gv->MPID == 0)
                  saveseis_glob(gv->SEISMO_FULLDATA, acq->recpos, acq->srcpos, ishot, gv->NS, 6, gv);

              break;
          default:
              break;

        }
    }
}
