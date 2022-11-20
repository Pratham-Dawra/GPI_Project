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

/* $Id: update_s_elastic_PML.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void update_s_visc_TTI_PML ( int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                            float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy,
                            float **sxx, float **syy, float **sxy,
                              float ***pr, float ***pp, float ***pq,
                            float ** pc11u, float **pc33u, float **pc13u, float ** pc55u, float ** pc15u, float ** pc35u,
                           float ** pc55ipjpu, float ** pc15ipjpu,float ** pc35ipjpu,
                           float *** pc11d, float ***pc33d, float ***pc13d, float *** pc55d,
                           float *** pc15d, float *** pc35d,
                           float *** pc55ipjpd, float *** pc15ipjpd,float *** pc35ipjpd,
                             float *bip, float *cip, float *hc,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                            float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx )
{
	int i,j, h1;
	extern int FW;

	

	


	/* interior */
	/*for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}*/

    /* left boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {

			cpml_update_s_x ( i,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
        }
    }

    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

		}
	}

    /* right boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            
            h1 = ( i-nx2+2*FW );

            cpml_update_s_x ( h1,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
                           b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
        }
    }

    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
  
            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

        }
    }

    /* top boundary */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            cpml_update_s_y ( i,j,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

        }
    }
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

        }
    }

    /* bottom boundary */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            h1 = ( j-ny2+2*FW );
            cpml_update_s_y ( i,h1,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
            
        }
    }
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

        }
    }

    /* corners */

    /*left-top*/
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            cpml_update_s_x ( i,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
                           b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

            cpml_update_s_y ( i,j,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
        }
    }
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

        }
    }

    /*left-bottom*/
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {

            cpml_update_s_x ( i,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
                           b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

            h1 = ( j-ny2+2*FW );

            cpml_update_s_y ( i,h1,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
           
        }
    }
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            
            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

            
        }
    }

    /* right-top */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            h1 = ( i-nx2+2*FW );
            cpml_update_s_x ( h1,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
                           b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

            cpml_update_s_y ( i,j,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
        }
    }
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

            
        }
    }

    /* right-bottom */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {

            h1 = ( i-nx2+2*FW );
            cpml_update_s_x ( h1,j,&pvxx[j][i],&pvyx[j][i],K_x,a_x,
                           b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

            h1 = ( j-ny2+2*FW );
            cpml_update_s_y ( i,h1,&pvxy[j][i],&pvyy[j][i],K_y,a_y,
                           b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

        }
    }
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            wavefield_update_s_visc_TTI ( i,j,pvxx,pvyx,pvxy,pvyy,sxy,sxx,syy,pr, pp, pq,
                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                         bip,  cip);

        }
    }


    
}
