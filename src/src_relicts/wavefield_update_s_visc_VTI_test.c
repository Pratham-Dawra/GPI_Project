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
                              float *bip, float *cip)
{
    extern float DT;
    extern int L;
    
    int l;
	float  dthalbe;
    float rv[L+1], pv[L+1], qv[L+1];
    float c55ipjpu, c11u, c13u, c33u;
    float c55ipjpd[L+1], c11d[L+1], c13d[L+1], c33d[L+1];

	
	float sumr=0.0, sump=0.0, sumq=0.0;
	/* computing sums of the old memory variables */
	
	dthalbe = DT/2.0;
    
	
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
        rv[l]=r[j][i][l];
        pv[l]=p[j][i][l];
        qv[l]=q[j][i][l];

		sumr+=rv[l];
		sump+=pv[l];
		sumq+=qv[l];
        
        c11d[l]=pc11d[j][i][l];
        c13d[l]=pc13d[j][i][l];
        c33d[l]=pc33d[j][i][l];
        c55ipjpd[l]=pc55ipjpd[j][i][l];

	}

    c55ipjpu=pc55ipjpu[j][i];
    c11u=pc11u[j][i];
    c13u=pc13u[j][i];
    c33u=pc33u[j][i];
    

    

    

	/* updating components of the stress tensor, partially */
	sxy[j][i] += ( c55ipjpu* ( vxy+vyx ) ) + ( dthalbe*sumr );
	sxx[j][i] += (c11u*vxx) + (c13u*vyy)+ ( dthalbe*sump );
	syy[j][i] += (c13u*vxx) + (c33u*vyy)+ ( dthalbe*sumq );



	/* now updating the memory-variables and sum them up*/
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
		rv[l] = bip[l]* ( rv[l]*cip[l]- ( c55ipjpd[l]* ( vxy+vyx ) ) );
		pv[l] = bip[l]* ( pv[l]*cip[l]- ( c11d[l]*vxx ) -  ( c13d[l]*vyy ) );
		qv[l] = bip[l]* ( qv[l]*cip[l]- ( c13d[l]*vxx ) -  ( c33d[l]*vyy ) );
		sumr += rv[l];
		sump += pv[l];
		sumq += qv[l];
        r[j][i][l]=rv[l];
        p[j][i][l]=pv[l];
        q[j][i][l]=qv[l];
	}


	/* and now the components of the stress tensor are
	   completely updated */
	sxy[j][i]+= ( dthalbe*sumr );
	sxx[j][i]+= ( dthalbe*sump );
	syy[j][i]+= ( dthalbe*sumq );
    
    

}
