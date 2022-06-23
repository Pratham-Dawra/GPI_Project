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

#include "fd.h"

void readmod_elastic_vti (float  **  rho, float **  pc11, float **  pc33, float **  pc13,float **  pc55 ){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, c11, c33, c13, c55;
    float vp, vs, epsilon, delta;
	int i, j, ii, jj;
    FILE *fp_epsilon, *fp_delta, *fp_vp, *fp_vs, *fp_rho;
    char filename[STRING_SIZE];

    
    fprintf(FP,"\n...reading elastic VTI Thomsen parameters from model-files...\n");

 fprintf(FP,"\n...reading viscoelastic TTI Thomsen parameters from model-files...\n");
 fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
 sprintf(filename,"%s.vp",MFILE);
 fp_vp=fopen(filename,"r");
 if ((fp_vp==NULL) && (MYID==0)) declare_error(" Could not open model file VP ! ");

  fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
  sprintf(filename,"%s.vs",MFILE);
  fp_vs=fopen(filename,"r");
  if ((fp_vs==NULL) && (MYID==0)) declare_error(" Could not open model file VS ! ");
  
  fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
  sprintf(filename,"%s.rho",MFILE);
  fp_rho=fopen(filename,"r");
  if ((fp_rho==NULL) && (MYID==0)) declare_error(" Could not open model file for densities ! ");


  fprintf(FP,"\t Epsilon:\n\t %s.epsilon\n\n",MFILE);
  sprintf(filename,"%s.epsilon",MFILE);
  fp_epsilon=fopen(filename,"r");
  if ((fp_epsilon==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Epsilon ! ");

  fprintf(FP,"\t Delta:\n\t %s.delta\n\n",MFILE);
  sprintf(filename,"%s.delta",MFILE);
  fp_delta=fopen(filename,"r");
  if ((fp_delta==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Delta ! ");



 

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
            for (j=1;j<=NYG;j++){

            fread(&vp, sizeof(float), 1, fp_vp);
            fread(&vs, sizeof(float), 1, fp_vs);
            fread(&rhov, sizeof(float), 1, fp_rho);
            fread(&epsilon, sizeof(float), 1, fp_epsilon);
            fread(&delta, sizeof(float), 1, fp_delta);

            c33=rhov*vp*vp;
            c55=rhov*vs*vs;
            c11=c33*(2.0*epsilon+1.0);
            c13=sqrt((2.0*delta*c33*(c33-c55))+((c33-c55)*(c33-c55)))-c55;




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
	



    fclose(fp_vp);
    fclose(fp_vs);
    fclose(fp_rho);
    fclose(fp_epsilon);
    fclose(fp_delta);


}




