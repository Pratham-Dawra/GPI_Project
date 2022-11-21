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

void readmod_elastic_vti (float  **  rho, float **  pc11, float **  pc33, float **  pc13,float **  pc55, GlobVar *gv ) {

	/* local variables */
	float rhov, c11, c33, c13, c55;
    float vp, vs, epsilon, delta;
	int i, j, ii, jj;
    FILE *fp_epsilon, *fp_delta, *fp_vp, *fp_vs, *fp_rho;
    char filename[STRING_SIZE+16];

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

    fprintf(gv->FP,"\n...reading elastic VTI Thomsen parameters from model-files...\n");

 fprintf(gv->FP,"\n...reading viscoelastic TTI Thomsen parameters from model-files...\n");
 fprintf(gv->FP,"\t Vp:\n\t %s.vp\n\n",gv->MFILE);
 sprintf(filename,"%s.vp",gv->MFILE);
 fp_vp=fopen(filename,"r");
 if ((fp_vp==NULL) && (MYID==0)) declare_error(" Could not open model file VP ! ");

  fprintf(gv->FP,"\t Vs:\n\t %s.vs\n\n",gv->MFILE);
  sprintf(filename,"%s.vs",gv->MFILE);
  fp_vs=fopen(filename,"r");
  if ((fp_vs==NULL) && (MYID==0)) declare_error(" Could not open model file VS ! ");
  
  fprintf(gv->FP,"\t Density:\n\t %s.rho\n\n",gv->MFILE);
  sprintf(filename,"%s.rho",gv->MFILE);
  fp_rho=fopen(filename,"r");
  if ((fp_rho==NULL) && (MYID==0)) declare_error(" Could not open model file for densities ! ");


  fprintf(gv->FP,"\t Epsilon:\n\t %s.epsilon\n\n",gv->MFILE);
  sprintf(filename,"%s.epsilon",gv->MFILE);
  fp_epsilon=fopen(filename,"r");
  if ((fp_epsilon==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Epsilon ! ");

  fprintf(gv->FP,"\t Delta:\n\t %s.delta\n\n",gv->MFILE);
  sprintf(filename,"%s.delta",gv->MFILE);
  fp_delta=fopen(filename,"r");
  if ((fp_delta==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Delta ! ");

	/* loop over global grid */
		for (i=1;i<=gv->NXG;i++){
            for (j=1;j<=gv->NYG;j++){

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
				if ((gv->POS[1]==((i-1)/gv->NX)) && 
				    (gv->POS[2]==((j-1)/gv->NY))){
					ii=i-gv->POS[1]*gv->NX;
					jj=j-gv->POS[2]*gv->NY;

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




