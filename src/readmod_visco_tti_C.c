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

void readmod_visco_tti (float  **  rho, float **  pc11, float **  pc33, float **  pc13, float **  pc55,
                        float **  pc15, float **  pc35,
                        float **  ptau11, float **  ptau33, float **  ptau13, float **  ptau55,
                        float ** ptau15, float ** ptau35, float *  eta){

    extern float DT, *FL, TS, TAU;
    extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];
    extern int WRITE_MODELFILES;
    extern FILE *FP;

		
	/* local variables */
    float rhov, c11, c33, c13, c55, tau11, tau33, tau13, tau55, tau15, tau35;
    float *pts, sumc11, sumc13, sumc33, sumc55, ws, fc;
    float l1, l2, l12, l22, l14, l24, l13, l23;
    float a1, a3, a4, a5, a6, t;
    float c11t, c33t, c55t, c13t, c15t, c35t;

	int i, j, l, ii, jj;
    FILE *fp_c11, *fp_c33, *fp_c13, *fp_c55, *fp_rho;
    FILE *fp_t11, *fp_t33, *fp_t13, *fp_t55, *fp_theta;

	char filename[STRING_SIZE];


    
	   fprintf(FP,"\n...reading VTI viscoelastic constants and density from model-files...\n");
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

        fprintf(FP,"\t TAU11:\n\t %s.tau11\n\n",MFILE);
        sprintf(filename,"%s.tau11",MFILE);
        fp_t11=fopen(filename,"r");
        if ((fp_t11==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants TAU11 ! ");
        
        fprintf(FP,"\t TAU33:\n\t %s.tau33\n\n",MFILE);
        sprintf(filename,"%s.tau33",MFILE);
        fp_t33=fopen(filename,"r");
        if ((fp_t33==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants TAU33 ! ");
        
        fprintf(FP,"\t TAU13:\n\t %s.tau13\n\n",MFILE);
        sprintf(filename,"%s.tau13",MFILE);
        fp_t13=fopen(filename,"r");
        if ((fp_t13==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants TAU13 ! ");
        
        fprintf(FP,"\t TAU55:\n\t %s.tau55\n\n",MFILE);
        sprintf(filename,"%s.tau55",MFILE);
        fp_t55=fopen(filename,"r");
        if ((fp_t55==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants TAU55 ! ");

        fprintf(FP,"\t Theta:\n\t %s.theta\n\n",MFILE);
        sprintf(filename,"%s.theta",MFILE);
        fp_theta=fopen(filename,"r");
        if ((fp_theta==NULL) && (MYID==0)) declare_error(" Could not open model file with elastic constants Theta ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if ((fp_rho==NULL) && (MYID==0)) declare_error(" Could not open model file for densities ! ");

	   
    /* vector for maxwellbodies */
    pts=vector(1,L);
    for (l=1;l<=L;l++) {
        pts[l]=1.0/(2.0*PI*FL[l]);
        eta[l]=DT/pts[l];
    }
   
    fc=1.0/TS;
   if (MYID==0){
        fprintf(FP," Message from readmod_visco_vti:\n");
        fprintf(FP," Center source frequency of %5.2f Hz applied for calculation of relaxed moduli ! \n",fc);
    }
    
    ws=2.0*PI*fc;
    


	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&c11, sizeof(float), 1, fp_c11);
			fread(&c33, sizeof(float), 1, fp_c33);
			fread(&c13, sizeof(float), 1, fp_c13);
			fread(&c55, sizeof(float), 1, fp_c55);
			fread(&rhov, sizeof(float), 1, fp_rho);
            fread(&tau11, sizeof(float), 1, fp_t11);
            fread(&tau33, sizeof(float), 1, fp_t33);
            fread(&tau13, sizeof(float), 1, fp_t13);
            fread(&tau55, sizeof(float), 1, fp_t55);
            fread(&t, sizeof(float), 1, fp_theta);

				
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
                if ((isnan(tau11)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in TAU11-Model !");}
                if ((isnan(tau13)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in TAU13-Model !");}
                if ((isnan(tau33)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in TAU33-Model !");}
                if ((isnan(tau55)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in TAU55-Model !");}
                if ((isnan(t)) && (MYID==0)) {
                                declare_error(" Found NaN-Values in Theta-Model !");}


                sumc11=0.0;  sumc33=0.0; sumc55=0.0;
                for (l=1;l<=L;l++){
                    sumc11=sumc11+((ws*ws*pts[l]*pts[l]*tau11)/(1.0+ws*ws*pts[l]*pts[l]));
                    sumc33=sumc33+((ws*ws*pts[l]*pts[l]*tau33)/(1.0+ws*ws*pts[l]*pts[l]));
                    sumc55=sumc55+((ws*ws*pts[l]*pts[l]*tau55)/(1.0+ws*ws*pts[l]*pts[l]));
                }

                /* relaxed moduli*/
                c11=c11/(1.0+sumc11);c33=c33/(1.0+sumc33);c55=c55/(1.0+sumc55);



                /* isotropic attenuation*/
                tau13=tau15=tau35=(c11*L*tau11-2.0*c55*L*tau55)/(c11-2.0*c55);

                sumc13=0.0;
                for (l=1;l<=L;l++){
                    sumc13=sumc13+((ws*ws*pts[l]*pts[l]*tau13)/(1.0+ws*ws*pts[l]*pts[l]));
                }
                c13=c13/(1.0+sumc13);
                
                
                /* Bond transformation (Oh et al, 2020, GJI, doi: 10.1093/gji/ggaa295 */
                t=t*PI/180.0;
                l1=cos(t); l2=sin(t); l12=l1*l1; l22=l2*l2; l14=l12*l12; l24=l22*l22; l13=l1*l12; l23=l2*l22;
               
                
                a1=2.0*c13+4.0*c55;
                a3=c11+c33-4.0*c55;
                a4=c11+c33-2.0*c13;
                a5=c13-c11+2.0*c55;
                a6=c13-c33+2.0*c55;
                
                c11t=c11*l14+c33*l24+a1*l12*l22;
                c33t=c11*l24+c33*l14+a1*l12*l22;
                c13t=a3*l12*l22+c13*(l14+l24);
                c55t=a4*l12*l22+c55*(l12-l22)*(l12-l22);
                c15t=a5*l13*l2-a6*l1*l23;
                c35t=a5*l23*l1-a6*l2*l13;

                
                

			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                    
                    pc11[jj][ii]=c11t;
                    rho[jj][ii]=rhov;
                    pc33[jj][ii]=c33t;
                    pc13[jj][ii]=c13t;
                    pc55[jj][ii]=c55t;
                    pc15[jj][ii]=c15t;
                    pc35[jj][ii]=c35t;
                    
                    ptau11[jj][ii]=tau11;
                    ptau33[jj][ii]=tau33;
                    ptau13[jj][ii]=tau13;
                    ptau55[jj][ii]=tau55;
                    ptau15[jj][ii]=tau15;
                    ptau35[jj][ii]=tau35;

				}
			}
		}
	




	fclose(fp_c11);
	fclose(fp_c13);
	fclose(fp_c33);
	fclose(fp_c55);
	fclose(fp_rho);
    fclose(fp_t11);
    fclose(fp_t13);
    fclose(fp_t33);
    fclose(fp_t55);
    
    
    
    
	
	/* each PE writes his model to disk */
	   
	   
/*	sprintf(filename,"%s.sofi2D.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	*/

}




