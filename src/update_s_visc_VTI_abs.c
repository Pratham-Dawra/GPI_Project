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
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* $Id: update_s_elastic_abs.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the absorbing frame (ABS=2 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void update_s_visc_vti_abs ( int *gx, int *gy, 
                              float **vx, float **vy, float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                             float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                             float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                             float *bip, float *cip, float ** absorb_coeff, float *hc, GlobVar *gv )
{
	
	int i,j,fdoh;
	float  vxx, vyy, vxy, vyx;


	fdoh=gv->FDORDER/2;

	/*Pointer array to the locations of the fd-operator functions*/
	void ( *FD_op_s[7] ) ();
	FD_op_s[1] = &operator_s_fd2;
	FD_op_s[2] = &operator_s_fd4;
	FD_op_s[3] = &operator_s_fd6;
	FD_op_s[4] = &operator_s_fd8;
	FD_op_s[5] = &operator_s_fd10;
	FD_op_s[6] = &operator_s_fd12;

	


	/* left boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			 wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
                                       bip,  cip, gv);

			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

 			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );

		}
	}

	/* top boundary */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* bottom boundary */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* corners */

	/*left-top*/
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );

		}
	}

	/*left-bottom*/
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right-top */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right-bottom */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc,gv );

			wavefield_update_s_visc_VTI ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pr, pp, pq,
	                    pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
	                                       bip,  cip, gv);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}


}
