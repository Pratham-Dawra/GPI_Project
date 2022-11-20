/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
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
/* -------------------------------------------------------------
 *  Precalculate derivates of particle velocities
 * -------------------------------------------------------------*/

#include "fd.h"

void v_derivatives(float **  vx, float **   vy, float **  pvxx, float **   pvyy, float **  pvyx, float **   pvxy, float *hc ){

	extern int NX, NY, FDORDER;
	int i, j, fdoh;
    float  vxx, vyy, vxy, vyx;

    
    
    fdoh=FDORDER/2;
    /*Pointer array to the locations of the fd-operator functions*/
    void ( *FD_op_s[7] ) ();
    FD_op_s[1] = &operator_s_fd2;
    FD_op_s[2] = &operator_s_fd4;
    FD_op_s[3] = &operator_s_fd6;
    FD_op_s[4] = &operator_s_fd8;
    FD_op_s[5] = &operator_s_fd10;
    FD_op_s[6] = &operator_s_fd12;
    
    for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            pvxx[j][i]=vxx;
            pvyy[j][i]=vyy;
            pvyx[j][i]=vyx;
            pvxy[j][i]=vxy;
		}
	}
}
