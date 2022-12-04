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
#include "logging.h"

void update_s_elastic_PML_4 ( int nx2, int ny2, int * gx, int * gy, int nt,
			      float **  vx, float **   vy, float **   sxx, float **   syy,
			      float **   sxy, float ** pi, float ** u, float ** uipjp,
			      float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
			      float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
			      float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx,float ** vxx_1,
			      float ** vxx_2,float ** vxx_3,float ** vxx_4,float ** vyy_1,float ** vyy_2,float ** vyy_3,
			      float ** vyy_4,float ** vxy_1,float ** vxy_2,float ** vxy_3,float ** vxy_4,float ** vyx_1,
			      float ** vyx_2,float ** vyx_3,float ** vyx_4, GlobVar *gv)
{
  int i,j, h1;
  float  vxx, vyy, vxy, vyx;
  double time1=0.0, time2=0.0;
  
  if ( ( gv->MPID==0 ) && ( ( nt+ ( gv->OUTNTIMESTEPINFO-1 ) ) %gv->OUTNTIMESTEPINFO ) ==0 ) {
    time1=MPI_Wtime();
    log_debug("Updating stress components...\n");
  }
    
  /* left boundary */
  for ( j=gy[2]+1; j<=gy[3]; j++ ) {
    for ( i=gx[1]; i<=gx[2]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* right boundary */
  for ( j=gy[2]+1; j<=gy[3]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      h1 = ( i-nx2+2*gv->FW );
      cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* top boundary */
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[2]+1; i<=gx[3]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y,b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* bottom boundary */
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[2]+1; i<=gx[3]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      h1 = ( j-ny2+2*gv->FW );
      cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* corners */
  
  /*left-top*/
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[1]; i<=gx[2]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y,b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /*left-bottom*/
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[1]; i<=gx[2]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      h1 = ( j-ny2+2*gv->FW );
      cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* right-top */
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      h1 = ( i-nx2+2*gv->FW );
      cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y, b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  /* right-bottom */
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      h1 = ( i-nx2+2*gv->FW );
      cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );
      h1 = ( j-ny2+2*gv->FW );
      cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );
      wavefield_update_s_el_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, gv );
    }
  }
    
  if ( ( gv->MPID==0 ) && ( ( nt+ ( gv->OUTNTIMESTEPINFO-1 ) ) %gv->OUTNTIMESTEPINFO ) ==0 ) {
    time2=MPI_Wtime();
    log_debug("Finished updating stress components (real time: %4.3fs).\n",time2-time1);
  }
}
