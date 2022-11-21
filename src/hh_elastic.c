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
 *   Model homogeneous half space
 *   if variable "h" is decreased, a layer over half-space is gained
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u, GlobVar *gv){

	/* local variables */
	float muv, piv, Vp, Vs, Rhov;
	float y;
	int i, j, ii, jj;
	char modfile[STRING_SIZE+16];
	float ** pwavemod=NULL, ** swavemod=NULL;

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	/*-----------------material property definition -------------------------*/	

	/* parameters for layer 1 */
	const float vp1=0.0, vs1=1.0, rho1=1.0, h=100.0;

	/* parameters for layer 2 */
	const float vp2=3500.0, vs2=2000.0, rho2=2000.0;

	/*-----------------------------------------------------------------------*/

	if (gv->WRITE_MODELFILES==1) {
		pwavemod  =  matrix(0,gv->NY+1,0,gv->NX+1);
		swavemod  =  matrix(0,gv->NY+1,0,gv->NX+1);
	}

	/* loop over global grid */
	for (i=1;i<=gv->NXG;i++){
		for (j=1;j<=gv->NYG;j++){

			/* calculate coordinate in m */
			y=(float)j*gv->DH;

			/* two layer case */
			if (y<=h){
				Vp=vp1; Vs=vs1; Rhov=rho1; }

			else{
				Vp=vp2; Vs=vs2; Rhov=rho2; }

			/* homogenous case */
			// 				vp=vp1; vs=vs1; rhov=rho1;

			muv=Vs*Vs*Rhov;
			piv=Vp*Vp*Rhov;

			/* only the PE which belongs to the current global gridpoint
				  is saving model parameters in his local arrays */
			if ((gv->POS[1]==((i-1)/gv->NX)) &&
					(gv->POS[2]==((j-1)/gv->NY))){
				ii=i-gv->POS[1]*gv->NX;
				jj=j-gv->POS[2]*gv->NY;

				u[jj][ii]=muv;
				rho[jj][ii]=Rhov;
				pi[jj][ii]=piv;
				if (gv->WRITE_MODELFILES==1)
				{
					pwavemod[jj][ii]=Vp;
					swavemod[jj][ii]=Vs;
				}
			}
		}
	}

	/* each PE writes his model to disk */

	/* only the density model is written to file */
	if (gv->WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI2D.rho",gv->MFILE);
		writemod(modfile,rho,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);
	}

	/* all models are written to file */
	if (gv->WRITE_MODELFILES==1) {
		sprintf(modfile,"%s.SOFI2D.u",gv->MFILE);
		writemod(modfile,u,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);

		sprintf(modfile,"%s.SOFI2D.pi",gv->MFILE);
		writemod(modfile,pi,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);

		sprintf(modfile,"%s.SOFI2D.vp",gv->MFILE);
		writemod(modfile,pwavemod,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);

		sprintf(modfile,"%s.SOFI2D.vs",gv->MFILE);
		writemod(modfile,swavemod,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);

		sprintf(modfile,"%s.SOFI2D.rho",gv->MFILE);
		writemod(modfile,rho,3,gv);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3,gv);
	}

	if (gv->WRITE_MODELFILES==1) {
		free_matrix(pwavemod,0,gv->NY+1,0,gv->NX+1);
		free_matrix(swavemod,0,gv->NY+1,0,gv->NX+1);
	}
}
