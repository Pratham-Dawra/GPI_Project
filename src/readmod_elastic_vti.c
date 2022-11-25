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
#include "logging.h"

void readmod_elastic_vti(float **rho, float **pc11, float **pc33, float **pc13,float **pc55, GlobVar *gv) 
{

  /* local variables */
  float rhov, c11, c33, c13, c55;
  float vp, vs, epsilon, delta;
  int ii, jj;
  char filename[STRING_SIZE+16];

  log_infoc(0, "Reading P-wave velocity model: %s.vp\n", gv->MFILE);
  sprintf(filename,"%s.vp",gv->MFILE);
  FILE *fp_vp=fopen(filename,"r");
  if (!fp_vp) log_fatal("Could not open P-wave velocity model %s.\n", filename);

  log_infoc(0, "Reading S-wave velocity model: %s.vs\n",gv->MFILE);
  sprintf(filename,"%s.vs",gv->MFILE);
  FILE *fp_vs=fopen(filename,"r");
  if (!fp_vs) log_fatal("Could not open S-wave velocity model %s.\n", filename);

  log_infoc(0, "Reading density model: %s.rho\n",gv->MFILE);
  sprintf(filename,"%s.rho",gv->MFILE);
  FILE *fp_rho=fopen(filename,"r");
  if (!fp_rho) log_fatal("Could not open density model %s\n", filename); 

  log_infoc(0, "Reading epsilon model: %s.epsilon\n",gv->MFILE);
  sprintf(filename,"%s.epsilon",gv->MFILE);
  FILE *fp_epsilon=fopen(filename,"r");
  if (!fp_epsilon) log_fatal("Could not open epsilon model %s.\n", filename);

  log_infoc(0, "Reading delta model: %s.delta\n",gv->MFILE);
  sprintf(filename,"%s.delta",gv->MFILE);
  FILE *fp_delta=fopen(filename,"r");
  if (!fp_delta) log_fatal("Could not open delta model %s.\n", filename);

  /* loop over global grid */
  for (int i=1;i<=gv->NXG;i++){
    for (int j=1;j<=gv->NYG;j++){
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

  return;
}




