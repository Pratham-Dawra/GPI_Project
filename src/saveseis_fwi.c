
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
 *   Write filtered seismograms and differences to file(s)
 * ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_fwi(int ishot, AcqVar *acq, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    float **fulldata = NULL;

    fulldata = matrix(1, gv->NTRG, 1, gv->NS);

    /* Write differences between measured and synthetic seismogramms (adjoint sources) to disk */
    if (vinv->WRITE_DIFF) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxdiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 11, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvydiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 12, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpdiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 14, gv);
        }
    }
    /* Write measured filtered seismogramms to disk */
    if (vinv->TIME_FILT && vinv->WRITE_FILTERED_DATA) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxdata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 21, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvydata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 22, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpdata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 24, gv);
        }
    }
    /* Write synthetic filtered seismogramms to disk */
    if (vinv->TIME_FILT && vinv->WRITE_FILTERED_DATA == 2) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxcalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 31, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvycalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 32, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpcalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 34, gv);
        }
    }

    free_matrix(fulldata, 1, gv->NTRG, 1, gv->NS);

}
