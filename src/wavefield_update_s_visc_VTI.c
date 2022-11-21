/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/* $Id: wavefield_update_s_visc.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the viscoelastic case*/

#include "fd.h"

void wavefield_update_s_visc_VTI ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy,
                              float **sxy, float **sxx, float ** syy, float ***p, float ***r,float ***q,
                              float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                              float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                              float *bip, float *cip, GlobVar *gv)
{
	int l;
	float  dthalbe;
	float sumr=0.0, sump=0.0, sumq=0.0;
	/* computing sums of the old memory variables */
	
	dthalbe = gv->DT/2.0;
	
	sumr=sump=sumq=0.0;
	for ( l=1; l<=gv->L; l++ ) {
		sumr+=r[j][i][l];
		sump+=p[j][i][l];
		sumq+=q[j][i][l];
	}


    
	/* updating components of the stress tensor, partially */
	sxy[j][i] += ( pc55ipjpu[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
	sxx[j][i] += (pc11u[j][i]*vxx) + (pc13u[j][i]*vyy)+ ( dthalbe*sump );
	syy[j][i] += (pc13u[j][i]*vxx) + (pc33u[j][i]*vyy)+ ( dthalbe*sumq );



	/* now updating the memory-variables and sum them up*/
	sumr=sump=sumq=0.0;
	for ( l=1; l<=gv->L; l++ ) {
		r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( pc55ipjpd[j][i][l]* ( vxy+vyx ) ) );
		p[j][i][l] = bip[l]* ( p[j][i][l]*cip[l]- ( pc11d[j][i][l]*vxx ) -  ( pc13d[j][i][l]*vyy ) );
		q[j][i][l] = bip[l]* ( q[j][i][l]*cip[l]- ( pc13d[j][i][l]*vxx ) -  ( pc33d[j][i][l]*vyy ) );
		sumr += r[j][i][l];
		sump += p[j][i][l];
		sumq += q[j][i][l];
	}


	/* and now the components of the stress tensor are
	   completely updated */
	sxy[j][i]+= ( dthalbe*sumr );
	sxx[j][i]+= ( dthalbe*sump );
	syy[j][i]+= ( dthalbe*sumq );

}
