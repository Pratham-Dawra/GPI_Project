/*------------------------------------------------------------------------
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
--------------------------------------------------------------------------*/

/* $Id: update_s_visc_interior.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at interior gridpoints (excluding boundarys) [gx2+1...gx3][gy2+1...gy3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_s_visc_interior ( int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                              float **vx, float **vy, float **sxx, float **syy,
                              float **sxy, float ***r, float *** p, float ***q,
                              float ** fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                              float *cjm, float ***d, float ***e, float ***dip,
                              float *hc, GlobVar *gv ) {

	int i,j,fdoh;
	float  vxx, vyy, vxy, vyx;
	double time1=0.0, time2=0.0;

    fdoh = gv->FDORDER/2;
    
    /*Pointer array to the locations of the fd-operator functions*/
    void ( *FD_op_s[7] ) ();
    FD_op_s[1] = &operator_s_fd2;
    FD_op_s[2] = &operator_s_fd4;
    FD_op_s[3] = &operator_s_fd6;
    FD_op_s[4] = &operator_s_fd8;
    FD_op_s[5] = &operator_s_fd10;
    FD_op_s[6] = &operator_s_fd12;
    

	if ( ( gv->MPID==0 ) && ( ( nt+ ( gv->OUTNTIMESTEPINFO-1 ) ) %gv->OUTNTIMESTEPINFO ) ==0 ) {
		time1=MPI_Wtime();
		log_debug("Updating stress components...\n");
	}


		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				
                FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc, gv );
                
                wavefield_update_s_visc ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                              fipjp,f,g,bip,bjm,cip,cjm,d,e,dip, gv );
				
			}
		}
		
	if ( ( gv->MPID==0 ) && ( ( nt+ ( gv->OUTNTIMESTEPINFO-1 ) ) %gv->OUTNTIMESTEPINFO ) ==0 ) {
		time2=MPI_Wtime();
		log_debug("Finished updating stess components (real time: %4.3fs).\n",time2-time1 );
	}
}
