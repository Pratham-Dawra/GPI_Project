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
 *   Filter seismograms in time domain with a Butterworth filter
 *   Lowpass or highpass filtering can be applied                                  
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"
#include "cseife.h"

void  timedomain_filt(float ** data, float vinv->F_LOW_PASS, int vinv->ORDER, int gv->NTRG, int gv->NS, int method){

	/* 
	data	: 	2-dimensional array containing seismograms (
	vinv->F_LOW_PASS	:	corner frequency in Hz
	vinv->ORDER	:	order of filter
	gv->NTRG	:	number of traces
	gv->NS	:	number of samples
	method	:	definition of filter
			1: lowpass filter
			2: highpass filter
	*/
	
	/* declaration of local variables */
	int itr, j;
	double *seismogram, T0;
	double *seismogram_hp, T0_hp;
	
	seismogram = dvector(1,gv->NS);
	
	seismogram_hp = dvector(1,gv->NS);
	
	T0=1.0/(double)vinv->F_LOW_PASS;
	if(vinv->F_HIGH_PASS)
		T0_hp=1.0/(double)vinv->F_HIGH_PASS;
	if(method==2)
		T0_hp=1.0/(double)vinv->F_LOW_PASS;
		
	if (method==1){    /* lowpass filter */
		for (itr=1;itr<=gv->NTRG;itr++){
			for (j=1;j<=gv->NS;j++){
				seismogram[j]=(double)data[itr][j];
			}
			
			seife_lpb(seismogram,gv->NS+1,gv->DT,T0,vinv->ORDER); /* gv->NS+1 because vector[0] is also allocated and otherwise seife_lpb do not filter the last sample */

			for (j=1;j<=gv->NS;j++){
				data[itr][j]=(float)seismogram[j];
			}
		}
	} /* end of itr<=gv->NTRG loop */

	if ((method==2)||(vinv->F_HIGH_PASS)){   /*highpass filter*/
		for (itr=1;itr<=gv->NTRG;itr++){
			for (j=1;j<=gv->NS;j++){
				seismogram_hp[j]=(double)data[itr][j];
			}
			
			seife_hpb(seismogram_hp,gv->NS+1,gv->DT,T0_hp,vinv->ORDER);
			
			for (j=1;j<=gv->NS;j++){
				data[itr][j]=(float)seismogram_hp[j];
			}
		}
	} /* end of itr<=gv->NTRG loop */
	
	free_dvector(seismogram,1,gv->NS);

	free_dvector(seismogram_hp,1,gv->NS);

} /* end of function */
