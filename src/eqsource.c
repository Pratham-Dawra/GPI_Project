/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
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
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Source excitation using Moment tensor for simulation of earthquakes (double-couple)
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void eqsource(int nt, AcqVar *acq, MemWavefield *mpw, GlobVar *gv)
{
	int i, j, l;
	float amp,dip,scale_amp;
	float m11,m22,m12;
		
	/* adding source wavelet to stress components
           (moment tensor source) at source points */
	for (l=1;l<=acq->nsrc_loc;l++) {
		i=(int)acq->srcpos_loc[1][l];
		j=(int)acq->srcpos_loc[2][l];

		dip = acq->srcpos_loc[7][l]*PI/180.0;
		scale_amp = gv->DT / (gv->DH * gv->DH);

		amp = acq->srcpos_loc[6][l] * scale_amp * acq->signals[l][nt];

		m11 = - sin(2*dip);
		m12 = - cos(2*dip);
		m22 =   sin(2*dip);		
		
		mpw->psxx[j][i] += amp * m11;
		mpw->psxy[j][i] += 0.3 * amp * m12;
		mpw->psxy[j-1][i] += 0.3 * amp * m12;
		mpw->psxy[j][i-1] += 0.3 * amp * m12;
		mpw->psyy[j][i] += amp * m22;
	}
}

