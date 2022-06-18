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
 *
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_visco_vti(float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55,
	float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, DH, TAU;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern int WRITE_MODELFILES;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float c11, c33, c55, c13, Rho, tau13;
	float *pts, sumc11, sumc13, sumc33, sumc55, ws;
	int i, j, l, ii, jj;
	char modfile[STRING_SIZE];	

	/*-----------------material property definition -------------------------*/	

	
/* anisotropic case, values from
                 Jones, Wang, 1981, Geophysics, 46, 3, 288-297*/
                
      const float C11=34.3e9, C33=22.7e9, C55=5.4e9, C13=10.7e9, RHO=2000.0;
      const float tau11=TAU, tau33=TAU, tau55=TAU;

      const int jh=50;

	/*-----------------------------------------------------------------------*/



	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];





	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){

			if (j<jh){ c11=0.0; c33=0.0; c13=0.0; c55=0.0; Rho=1.0;}
			else { c11=C11; c33=C33; c13=C13; c55=C55; Rho=RHO;}
	

			sumc11=0.0;  sumc33=0.0; sumc55=0.0;
			for (l=1;l<=L;l++){
				sumc11=sumc11+((ws*ws*pts[l]*pts[l]*tau11)/(1.0+ws*ws*pts[l]*pts[l]));
				sumc33=sumc33+((ws*ws*pts[l]*pts[l]*tau33)/(1.0+ws*ws*pts[l]*pts[l]));
				sumc55=sumc55+((ws*ws*pts[l]*pts[l]*tau55)/(1.0+ws*ws*pts[l]*pts[l]));
			}

			/* relaxed moduli*/
			c11=c11/(1.0+sumc11);c33=c33/(1.0+sumc33);c55=c55/(1.0+sumc55);



			/* isotropic attenuation*/
			tau13=(c11*L*tau11-2.0*c55*L*tau55)/(c11-2.0*c55);
			sumc13=0.0;
			for (l=1;l<=L;l++){
				sumc13=sumc13+((ws*ws*pts[l]*pts[l]*tau13)/(1.0+ws*ws*pts[l]*pts[l]));
			}
			c13=c13/(1.0+sumc13);


			

			/* only the PE which belongs to the current global gridpoint
				  is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) &&
					(POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				pc11[jj][ii]=c11;
				rho[jj][ii]=Rho;
				pc33[jj][ii]=c33;
				pc13[jj][ii]=c13;
				pc55[jj][ii]=c55;

				ptau11[jj][ii]=tau11;
				ptau33[jj][ii]=tau33;
				ptau13[jj][ii]=tau13;
				ptau55[jj][ii]=tau55;


			}
		}
	}




	/* each PE writes his model to disk */

		/* all models are written to file */
	if (WRITE_MODELFILES==1) {
		sprintf(modfile,"%s.SOFI2D.c11",MFILE);
		writemod(modfile,pc11,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.c33",MFILE);
		writemod(modfile,pc33,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.c13",MFILE);
		writemod(modfile,pc13,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.c55",MFILE);
		writemod(modfile,pc55,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);


		sprintf(modfile,"%s.SOFI2D.tau11",MFILE);
		writemod(modfile,ptau11,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.tau33",MFILE);
		writemod(modfile,ptau33,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.tau13",MFILE);
		writemod(modfile,ptau13,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.tau55",MFILE);
		writemod(modfile,ptau55,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	/* only density is written to file */
	if (WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	free_vector(pts,1,L);
}
