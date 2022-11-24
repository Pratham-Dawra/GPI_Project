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
/* ----------------------------------------------------------------------
 * Computation of local receiver coordinates
 * (within each subgrid)
 *	 
 * ----------------------------------------------------------------------*/
 
#include "fd.h"

int **splitrec(int **recpos, int *ntr_loc, int ntr, int *recswitch, GlobVar *gv) 
{
  int a,b,i=0,j,k;
  int **recpos_dummy = NULL, **recpos_local = NULL;

  recpos_dummy = imatrix(1,3,1,ntr);

  for (j=1;j<=ntr;j++) {
    recswitch[j] = 0;
    a=(recpos[1][j]-1)/gv->IENDX;
    b=(recpos[2][j]-1)/gv->IENDY;
    if ((gv->POS[1]==a)&&(gv->POS[2]==b)) {
      recswitch[j] = 1;
      i++;	/* number of receivers i of each process */
      recpos_dummy[1][i] = ((recpos[1][j]-1)%gv->IENDX)+1;
      recpos_dummy[2][i] = ((recpos[2][j]-1)%gv->IENDY)+1;
      recpos_dummy[3][i] = j;
    }
  }
   
  if (i>0) recpos_local = imatrix(1,3,1,i);
	
  for (k=1;k<=i;k++){
    recpos_local[1][k] = recpos_dummy[1][k];
    recpos_local[2][k] = recpos_dummy[2][k];
    recpos_local[3][k] = recpos_dummy[3][k];
  }
	
  free_imatrix(recpos_dummy,1,3,1,ntr);
  
  *ntr_loc=i;

  return recpos_local;
}
