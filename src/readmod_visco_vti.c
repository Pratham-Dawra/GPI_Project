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

void readmod_visco_vti (float  **  rho,
                            float **  pc11, float **  pc33, float **  pc13,float **  pc55,
                            float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55,
                        float *  eta, GlobVar *gv) {
    
	/* local variables */
    float rhov, c11, c33, c13, c55, tau11, tau33, tau13, tau55;
    float *pts, sumc11, sumc13, sumc33, sumc55, ws, fc;
    float vp, vs, qp, qs, epsilon, delta;
 
	int i, j, l, ii, jj;
    FILE *fp_epsilon, *fp_delta, *fp_vp, *fp_vs;
    FILE *fp_qp, *fp_qs, *fp_rho;

	char filename[STRING_SIZE+16];

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    
    fprintf(gv->FP,"\n...reading viscoelastic VTI Thomsen parameters from model-files...\n");
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

     fprintf(gv->FP,"\t Qp:\n\t %s.qp\n\n",gv->MFILE);
     sprintf(filename,"%s.qp",gv->MFILE);
     fp_qp=fopen(filename,"r");
     if ((fp_qp==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Qp ! ");

     fprintf(gv->FP,"\t Qs:\n\t %s.qs\n\n",gv->MFILE);
     sprintf(filename,"%s.qs",gv->MFILE);
     fp_qs=fopen(filename,"r");
     if ((fp_qs==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Qs ! ");

    /* vector for maxwellbodies */
    pts=vector(1,gv->L);
    for (l=1;l<=gv->L;l++) {
        pts[l]=1.0/(2.0*PI*gv->FL[l]);
        eta[l]=gv->DT/pts[l];
    }

    fc=1.0/gv->TS;
   if (MYID==0){
        fprintf(gv->FP," Message from readmod_visco_vti:\n");
        fprintf(gv->FP," Center source frequency of %5.2f Hz applied for calculation of relaxed moduli ! \n",fc);
    }

    ws=2.0*PI*fc;

	/* loop over global grid */
		for (i=1;i<=gv->NXG;i++){
			for (j=1;j<=gv->NYG;j++){
			
                fread(&vp, sizeof(float), 1, fp_vp);
                fread(&vs, sizeof(float), 1, fp_vs);
                fread(&qp, sizeof(float), 1, fp_qp);
                fread(&qs, sizeof(float), 1, fp_qs);
                fread(&rhov, sizeof(float), 1, fp_rho);
                fread(&epsilon, sizeof(float), 1, fp_epsilon);
                fread(&delta, sizeof(float), 1, fp_delta);

                c33=rhov*vp*vp;
                c55=rhov*vs*vs;
                c11=c33*(2.0*epsilon+1.0);
                c13=sqrt((2.0*delta*c33*(c33-c55))+((c33-c55)*(c33-c55)))-c55;
                tau11=tau33=2.0/(qp*gv->L); tau55=2.0/(qs*gv->L);
 
                sumc11=0.0;  sumc33=0.0; sumc55=0.0;
                for (l=1;l<=gv->L;l++){
                    sumc11=sumc11+((ws*ws*pts[l]*pts[l]*tau11)/(1.0+ws*ws*pts[l]*pts[l]));
                    sumc33=sumc33+((ws*ws*pts[l]*pts[l]*tau33)/(1.0+ws*ws*pts[l]*pts[l]));
                    sumc55=sumc55+((ws*ws*pts[l]*pts[l]*tau55)/(1.0+ws*ws*pts[l]*pts[l]));
                }

                /* relaxed moduli*/
                c11=c11/(1.0+sumc11);c33=c33/(1.0+sumc33);c55=c55/(1.0+sumc55);

                /* isotropic attenuation*/
                tau13=(c11*gv->L*tau11-2.0*c55*gv->L*tau55)/(c11-2.0*c55);

                sumc13=0.0;
                for (l=1;l<=gv->L;l++){
                    sumc13=sumc13+((ws*ws*pts[l]*pts[l]*tau13)/(1.0+ws*ws*pts[l]*pts[l]));
                }
                c13=c13/(1.0+sumc13);

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
                    ptau11[jj][ii]=tau11;
                    ptau33[jj][ii]=tau33;
                    ptau13[jj][ii]=tau13;
                    ptau55[jj][ii]=tau55;
				}
			}
		}

    fclose(fp_vp);
    fclose(fp_vs);
    fclose(fp_qp);
    fclose(fp_qs);
    fclose(fp_rho);
    fclose(fp_epsilon);
    fclose(fp_delta);

    free_vector(pts,1,gv->L);

}




