/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* $Id: update_s_elastic_interior.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at interior gridpoints (excluding boundarys) [gx2+1...gx3][gy2+1...gy3]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void update_s_elastic_VTI_interior ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                        float **  vx, float **   vy, float **   sxx, float **   syy,
                        float **   sxy, float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33, float *hc )
{


	int i,j;
	float c55ipjp, c33, c11, c13;
	float  vxx, vyy, vxy, vyx;
	float  dhi;
	extern float DT, DH;
	extern int MYID, FDORDER;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;


	dhi = 1.0/DH;

	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time1=MPI_Wtime();
		fprintf ( FP,"\n **Message from update_s_interior (printed by PE %d):\n",MYID );
		fprintf ( FP," Updating stress components ..." );
	}

	if (FDORDER>2) declare_error("FDORDER>2 not available !");
	
	switch ( FDORDER ) {

	case 2:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {


				/* spatial derivatives of the components of the velocities */
				/* using Holberg coefficients */
				vxx = hc[1]* ( vx[j][i]  -vx[j][i-1] ) *dhi;
				vyy = hc[1]* ( vy[j][i]  -vy[j-1][i] ) *dhi;
				vyx = hc[1]* ( vy[j][i+1]-vy[j][i] ) *dhi;
				vxy = hc[1]* ( vx[j+1][i]-vx[j][i] ) *dhi;

				/* updating components of the stress tensor, partially */
				c55ipjp=pc55ipjp[j][i]*DT;
				c11=pc11[j][i]*DT;
				c33=pc33[j][i]*DT;
				c13=pc13[j][i]*DT;

				sxy[j][i]+= ( c55ipjp* ( vxy+vyx ) );

				sxx[j][i]+= ( (c11* vxx)+ (c13*vyy) );
				syy[j][i]+= ( (c13* vxx)+ (c33*vyy) );

/*				sxy[j][i]+= ( fipjp* ( vxy+vyx ) );
				sxx[j][i]+= ( g* ( vxx+vyy ) )- ( 2.0*f*vyy );
				syy[j][i]+= ( g* ( vxx+vyy ) )- ( 2.0*f*vxx );*/

			}
		}
		break;



	} /* end of switch(FDORDER) */



	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time2=MPI_Wtime();
		fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
	}
}


