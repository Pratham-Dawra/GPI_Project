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
#include "logging.h"

void update_s_elastic_abs ( int * gx, int * gy, int nt,
                            float **  vx, float **   vy, float **   sxx, float **   syy,
                            float **   sxy, float ** pi, float ** u, float ** uipjp,
                            float ** absorb_coeff, GlobVar *gv )
{
  int i,j;
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
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* right boundary */
  for ( j=gy[2]+1; j<=gy[3]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* top boundary */
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[2]+1; i<=gx[3]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* bottom boundary */
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[2]+1; i<=gx[3]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* corners */

  /*left-top*/
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[1]; i<=gx[2]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /*left-bottom*/
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[1]; i<=gx[2]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* right-top */
  for ( j=gy[1]; j<=gy[2]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  /* right-bottom */
  for ( j=gy[3]+1; j<=gy[4]; j++ ) {
    for ( i=gx[3]+1; i<=gx[4]; i++ ) {
      gv->FDOP_S( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy );
      wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp, gv );
      abs_update_s ( i,j,sxx,sxy,syy,absorb_coeff );
    }
  }

  if ( ( gv->MPID==0 ) && ( ( nt+ ( gv->OUTNTIMESTEPINFO-1 ) ) %gv->OUTNTIMESTEPINFO ) ==0 ) {
    time2=MPI_Wtime();
    log_debug("Finished updating stress components (real time: %4.3fs).\n",time2-time1);
  }
}
