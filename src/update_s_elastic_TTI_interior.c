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

void update_s_elastic_TTI_interior ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                                    float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                                    float **   sxx, float **   syy,
                        float **   sxy, float ** pc11, float ** pc55ipjp, float ** pc13, float ** pc33,
                                    float ** pc15, float ** pc35, float ** pc15ipjp, float ** pc35ipjp, float *hc )
{


	int i,j;
	extern int MYID;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;

  


    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time1=MPI_Wtime();
        fprintf ( FP,"\n **Message from update_s_elastic_tti_interior (printed by PE %d):\n",MYID );
        fprintf ( FP," Updating stress components ..." );
    }
    

    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {
            wavefield_update_s_el_tti (i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pc11,pc55ipjp,pc13,pc33, pc15, pc35, pc15ipjp, pc35ipjp);
        }
    }
    
    
    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time2=MPI_Wtime();
        fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
    }
}


