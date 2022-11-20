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

void update_s_visc_tti_abs ( int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                            float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                            float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                            float ** pc11u, float **pc33u, float **pc13u, float ** pc55u, float ** pc15u, float ** pc35u,
                           float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                           float *** pc11d, float ***pc33d, float ***pc13d, float *** pc55d,
                           float *** pc15d, float *** pc35d,
                           float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                             float *bip, float *cip, float ** absorb_coeff, float *hc )
{
	
	int i,j;
	


	
	


	/* left boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
            
            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
            
            abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );

		}
	}

	/* top boundary */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* bottom boundary */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* corners */

	/*left-top*/
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );

		}
	}

	/*left-bottom*/
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right-top */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}

	/* right-bottom */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);
			abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
		}
	}


}
