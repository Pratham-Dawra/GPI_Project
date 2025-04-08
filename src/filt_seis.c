
/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public LiceNTe as published by
 * the Free Software Foundation, version 2.0 of the LiceNTe only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public LiceNTe for more details.
 * 
 * You should have received a copy of the GNU General Public LiceNTe
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/liceNTes/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
* Filter seismograms in time domain with a Butterworth filter of the libseife
* library. Lowpass or highpass filtering can be applied                                 
* last update 2011, L. Rehor
*  ----------------------------------------------------------------------*/

#include "fd.h"
#include "util.h"
#include "cseife.h"

void filt_seis(float **data, int ntr, int NS, float finv, const GlobVar *gv)
{
    /*
     * data    :       2-dimensional array containing seismograms (
     * finv    :       corner frequency in Hz
     * order   :       order of filter
     * ntr     :       number of traces
     * NS      :       number of samples
     * method  :       definition of filter
     * 1: lowpass filter
     * 2: highpass filter
     */

    const int order = 4;
    // const int method = 1;
    float fc = finv;
    double *seismogram = dvector(0, NS-1);
    double T0 = 1.0 / fc;

    for (int itr = 1; itr <= ntr; itr++) {
        for (int j = 0; j <= NS-1; j++) {
	    seismogram[j] = (double)data[itr][j];
        }

        // if (method == 1) {      /* lowpass filter */
	seife_lpb(seismogram, NS, gv->DT, T0, order);
	// }

        // if (method == 2) {      /* highpass filter */
        //     seife_hpb(seismogram, NT, gv->DT, T0, order);
        // }

        for (int j = 0; j <= NS-1; j++) {
            data[itr][j] = (float)seismogram[j];
        }
    }

    free_dvector(seismogram, 0, NS-1);
}
