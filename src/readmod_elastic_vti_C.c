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
/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files
 *   This file contains function readmod, which has the purpose
 *   to read data from model-files for viscoelastic simulation
 *
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_elastic_vti (float  **  rho, float **  pc11, float **  pc33, float **  pc13,float **  pc55 ){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, c11, c33, c13, c55;
	int i, j, ii, jj;
	FILE *fp_c11, *fp_c33, *fp_c13, *fp_c55, *fp_rho;
	char filename[STRING_SIZE];




	   fprintf(FP,"\n...reading VTI elastic constants and density from model-files...\n");

	   fprintf(FP,"\t C11:\n\t %s.c11\n\n",MFILE);
	   sprintf(filename,"%s.c11",MFILE);
	   fp_c11=fopen(filename,"r");
	   if ((fp_c11==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants C11 ! ");

 		fprintf(FP,"\t C33:\n\t %s.c33\n\n",MFILE);
	   sprintf(filename,"%s.c33",MFILE);
	   fp_c33=fopen(filename,"r");
	   if ((fp_c33==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants C33 ! ");

		fprintf(FP,"\t C13:\n\t %s.c13\n\n",MFILE);
	   sprintf(filename,"%s.c13",MFILE);
	   fp_c13=fopen(filename,"r");
	   if ((fp_c13==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants C13 ! ");

		fprintf(FP,"\t C55:\n\t %s.c55\n\n",MFILE);
	   sprintf(filename,"%s.c55",MFILE);
	   fp_c55=fopen(filename,"r");
	   if ((fp_c55==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants C55 ! ");


	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if ((fp_rho==NULL) && (MYID==0)) declare_error(" Could not open model file for densities ! ");

	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&c11, sizeof(float), 1, fp_c11);
			fread(&c33, sizeof(float), 1, fp_c33);
			fread(&c13, sizeof(float), 1, fp_c13);
			fread(&c55, sizeof(float), 1, fp_c55);
			fread(&rhov, sizeof(float), 1, fp_rho);
				
		if ((isnan(c11)) && (MYID==0)) {
           	declare_error(" Found NaN-Values in C11-Model !");}
        if ((isnan(c13)) && (MYID==0)) {
           	declare_error(" Found NaN-Values in C13-Model !");}
        if ((isnan(c33)) && (MYID==0)) {
           	declare_error(" Found NaN-Values in C33-Model !");}
        if ((isnan(c55)) && (MYID==0)) {
           	declare_error(" Found NaN-Values in C55-Model !");}
		if ((isnan(rhov)) && (MYID==0)) {
                declare_error(" Found NaN-Values in Rho-Model !");}



			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					pc11[jj][ii]=c11;
					pc13[jj][ii]=c13;
					pc33[jj][ii]=c33;
					pc55[jj][ii]=c55;
					rho[jj][ii]=rhov;
				}
			}
		}
	




	fclose(fp_c11);
	fclose(fp_c13);
	fclose(fp_c33);
	fclose(fp_c55);
	fclose(fp_rho);
	
	
	/* each PE writes his model to disk */
	   
	   
/*	sprintf(filename,"%s.sofi2D.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	*/

}




