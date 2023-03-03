
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
 *   generate P-wave source at source nodes
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, AcqVar *acq, MemWavefield * mpw, GlobVar * gv)
{
    float amp = 0.0f;
    int i, j;

    /* adding source wavelet to stress components (explosive source) at source points */

    for (int l = 1; l <= acq->nsrc_loc; l++) {
        i = (int)acq->srcpos_loc[1][l];
        j = (int)acq->srcpos_loc[2][l];

        //amp=acq->signals[l][nt]; //unscaled explosive source
        amp = (acq->signals[l][nt]) / (gv->DH * gv->DH); //scaled explosive source, seismic Moment = 1 Nm

        if (nt == 1) {
            amp = acq->signals[l][nt + 1] / (2.0 * gv->DH * gv->DH);
        }
        if ((nt > 1) && (nt < gv->NT)) {
            amp = (acq->signals[l][nt + 1] - acq->signals[l][nt - 1]) / (2.0 * gv->DH * gv->DH);
        }
        if (nt == gv->NT) {
            amp = -acq->signals[l][nt - 1] / (2.0 * gv->DH * gv->DH);
        }

        mpw->psxx[j][i] += amp;
        mpw->psyy[j][i] += amp;
    }
}
