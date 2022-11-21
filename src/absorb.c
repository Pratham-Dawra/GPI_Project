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
 *   Create dissipative boundarie around the model grid
 *   The dissipative coefficients are stored in the 2-D array
 *   absorb_coeff. The interior of the model is weighted by the
 *   coefficient 1. In the absorbing frame the coefficients 
 *   are less than one. Coefficients are computed using 
 *   exponential damping (see Cerjan et al., 1985, 
 *   Geophysics, 50, 705-708)
 *
 ------------------------------------------------------------- */

#include "fd.h"

void absorb(float ** absorb_coeff, GlobVar *gv)
{

	/* extern variables */

	//extern float DAMPING;
	//extern int FREE_SURF, NX, NY, BOUNDARY, FW;
	//extern int NPROCX, NPROCY, MYID, POS[3];
	//extern FILE *FP;
	
	/* local variables */
	int i, j, ii, jj, xb, yb, xe, ye;
	float amp, a, *coeff;
	/*char modfile[STRING_SIZE];*/

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	if (MYID==0){
		fprintf(gv->FP,"\n **Message from absorb (printed by PE %d):\n",MYID);
		fprintf(gv->FP," Coefficients for absorbing frame are now calculated.\n");
		fprintf(gv->FP," Width of dissipative frame (grid points)= %i\n",gv->FW);
		fprintf(gv->FP," Percentage of exponential damping = %5.2f\n",gv->DAMPING);
	}

	amp=1.0-gv->DAMPING/100.0;   /* amplitude at the edge of the numerical grid */
	
	coeff=vector(1,gv->FW);
	a=sqrt(-log(amp)/((gv->FW-1)*(gv->FW-1)));
	
	for (i=1;i<=gv->FW;i++)
		coeff[i]=exp(-(a*a*(gv->FW-i)*(gv->FW-i)));
	
	if (MYID==0)
	{
		fprintf(gv->FP," Table of coefficients \n # \t coeff \n");
		/*printf(" gv->FW=%d \t a=%f amp=%f \n", gv->FW,a,amp); */
		for (i=1;i<=gv->FW;i++)
			fprintf(gv->FP," %d \t %5.3f \n", i, coeff[i]);
	}	
	

	/* initialize array of coefficients with one */
	for (j=1;j<=gv->NY;j++)
	for (i=1;i<=gv->NX;i++) absorb_coeff[j][i]=1.0;


	/* compute coefficients for left and right grid boundaries (x-direction) */
	if ((!gv->BOUNDARY) && (gv->POS[1]==0))
	{
		yb=1; ye=gv->NY; 
		for (i=1;i<=gv->FW;i++)
		{
			if ((gv->POS[2]==0) && (!(gv->FREE_SURF))) yb=i;
			if (gv->POS[2]==gv->NPROCY-1) ye=gv->NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][i]=coeff[i];
		}
	}
			
	if ((!gv->BOUNDARY) && (gv->POS[1]==gv->NPROCX-1))
	{
		yb=1; ye=gv->NY;
		for (i=1;i<=gv->FW;i++){
			ii=gv->NX-i+1;
			if ((gv->POS[2]==0) && (!(gv->FREE_SURF))) yb=i;
			if (gv->POS[2]==gv->NPROCY-1) ye=gv->NY-i+1;
			for (j=yb;j<=ye;j++)
				absorb_coeff[j][ii]=coeff[i];
		}
	}
	

	/* compute coefficients for top and bottom grid boundaries (y-direction) */

	if ((gv->POS[2]==0) && (!(gv->FREE_SURF)))
	{
		xb=1; xe=gv->NX;
		for (j=1;j<=gv->FW;j++)
		{
			if ((!gv->BOUNDARY) && (gv->POS[1]==0)) xb=j;
			if ((!gv->BOUNDARY) && (gv->POS[1]==gv->NPROCX-1)) xe=gv->NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[j][i]=coeff[j];
		}
	}

	if (gv->POS[2]==gv->NPROCY-1)
	{
		xb=1; xe=gv->NX;
		for (j=1;j<=gv->FW;j++)
		{
			jj=gv->NY-j+1;
			if ((!gv->BOUNDARY) && (gv->POS[1]==0)) xb=j;
			if ((!gv->BOUNDARY) && (gv->POS[1]==gv->NPROCX-1)) xe=gv->NX-j+1;
			for (i=xb;i<=xe;i++)
				absorb_coeff[jj][i]=coeff[j];
		}
	}


/*	sprintf(modfile,"absorb_coeff.bin");

	writemod(modfile,absorb_coeff,3); 

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3); 
*/

	free_vector(coeff,1,gv->FW);
}



