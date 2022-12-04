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

void v_derivatives(float **vx, float **vy, float **pvxx, float **pvyy, float **pvyx, float **pvxy, GlobVar *gv)
{
  float vxx, vyy, vxy, vyx;
  
  for (int j=1; j<=gv->NY; j++){
    for (int i=1; i<=gv->NX; i++){
      gv->FDOP_S(i,j,&vxx,&vyx,&vxy,&vyy,vx,vy);
      pvxx[j][i] = vxx;
      pvyy[j][i] = vyy;
      pvyx[j][i] = vyx;
      pvxy[j][i] = vxy;
    }
  }
}
