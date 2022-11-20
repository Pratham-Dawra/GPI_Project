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

void wavefield_update_s_visc_TTI ( int i, int j,float  **vxx, float **vyx,float **vxy,float **vyy,
                              float **sxy, float **sxx, float ** syy, float ***p, float ***r,
                                  float ***q,
                                  float ** pc11u, float **pc33u, float **pc13u, float ** pc55u, float ** pc15u, float ** pc35u,
                                 float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                                 float *** pc11d, float ***pc33d, float ***pc13d, float *** pc55d,
                                 float *** pc15d, float *** pc35d,
                                 float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                                  float *bip, float *cip)
{
	int l;
	float  dthalbe;
	extern float DT;
	extern int L;
	float sumr=0.0, sump=0.0, sumq=0.0;
    float vxxipjp, vyyipjp, vyxij, vxyij, vij, v;

    
    
    
	/* computing sums of the old memory variables */
	
	dthalbe = DT/2.0;
	
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
		sumr+=r[j][i][l];
		sump+=p[j][i][l];
		sumq+=q[j][i][l];
	}

        
    vxxipjp=0.25*(vxx[j][i]+vxx[j+1][i]+vxx[j][i+1]+vxx[j+1][i+1]);
    vyyipjp=0.25*(vyy[j][i]+vyy[j+1][i]+vyy[j][i+1]+vyy[j+1][i+1]);
    
    vyxij=0.25*(vyx[j][i]+vyx[j-1][i]+vyx[j][i-1]+vyx[j-1][i-1]);
    vxyij=0.25*(vxy[j][i]+vxy[j-1][i]+vxy[j][i-1]+vxy[j-1][i-1]);
    vij=vyxij+vxyij;

    v=vxy[j][i]+vyx[j][i];

    
	/* updating components of the stress tensor, partially */
	sxy[j][i] += ( pc55ipjpu[j][i]* v ) + ( pc15ipjpu[j][i]* vxxipjp ) + ( pc35ipjpu[j][i]* vyyipjp )+ ( dthalbe*sumr );
	sxx[j][i] += (pc11u[j][i]*vxx[j][i]) + (pc13u[j][i]*vyy[j][i])+ + (pc15u[j][i]*vij) + ( dthalbe*sump );
	syy[j][i] += (pc13u[j][i]*vxx[j][i]) + (pc33u[j][i]*vyy[j][i])+ (pc35u[j][i]*vij) + ( dthalbe*sumq );



	/* now updating the memory-variables and sum them up*/
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
		r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( pc55ipjpd[j][i][l]* v ) - ( pc15ipjpd[j][i][l]* vxxipjp) - ( pc35ipjpd[j][i][l]* vyyipjp) );
		p[j][i][l] = bip[l]* ( p[j][i][l]*cip[l]- ( pc11d[j][i][l]*vxx[j][i] ) -  ( pc13d[j][i][l]*vyy[j][i] ) - ( pc15d[j][i][l]*vij ) );
		q[j][i][l] = bip[l]* ( q[j][i][l]*cip[l]- ( pc13d[j][i][l]*vxx[j][i] ) -  ( pc33d[j][i][l]*vyy[j][i] ) - ( pc35d[j][i][l]*vij ));
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
