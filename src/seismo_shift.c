/*------------------------------------------------------------------------
 * Copyright (C) 2024 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 *------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * shift seismograms back towards time zero
 *----------------------------------------------------------------------*/

#include <assert.h>
#include "fd.h"
#include "seismo_shift.h"
#include "globvar_struct.h"
#include "macros.h"

void seismo_shift(st_seismogram *section,
		  int ntr,
		  float fc,
		  const GlobVar *gv)
{
  float ts = 0.0f;
  int nshift, k;

  if ((0 == gv->AUTO_SHIFT) && (0 == gv->NDTSHIFT)) {
    // nothing to do
    return;
  }

  float **seiscomp = NULL;

  for (int c = 1; c <= 6; c++) {
    switch (c) {
    case 1:
      seiscomp = section->vx;
      break;
    case 2:
      seiscomp = section->vy;
      break;
    case 3:
      seiscomp = section->vz;
      break;
    case 4:
      seiscomp = section->curl;
      break;
    case 5:
      seiscomp = section->div;
      break;
    case 6:
      seiscomp = section->p;
      break;
    default:
      // should never ever reach this point
      seiscomp = NULL;
      break;
    }

    // if component is not allocated, skip
    if (!seiscomp) continue;

    for (int n = 1; n <= ntr; ++n) {
      nshift = 0;
      if (gv->AUTO_SHIFT) {
        ts = 0.0f;
        if (gv->SOURCE_SHAPE == 1) {
          ts = 1.5 / fc;           // time shift of zero-phase Ricker wavelet
        }
        // convert time shift into sample shift (rounded)
        nshift += iround(ts / gv->DT);
      } 
      // add user-requested NDTSHIFT
      nshift += gv->NDTSHIFT;
      
      if (nshift == 0) continue;
      
      k = 0;
      // shift trace towards the front of the buffer
      for (int i = nshift; i <= gv->NS-1; i++) {
        seiscomp[n][k++] = seiscomp[n][i];
      }
      // fill remaining empty slots with zeros
      for (int i = k; i <= gv->NS-1; i++) {
        seiscomp[n][i] = 0.0f;
      }
    } // trace-loop
  } // component-loop

  return;
}
