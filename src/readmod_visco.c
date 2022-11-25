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

void readmod_visco(float **rho, float **pi, float **u, float **taus, float **taup, float *eta, GlobVar *gv) 
{
  /* local variables */
  float rhov, muv, piv, vp, vs, qp, qs;
  float *pts, ts, tp, sumu, sumpi, ws;
  int l, ii, jj;
  char filename[STRING_SIZE+16];

  /* vector for maxwell bodies */
  pts=vector(1,gv->L);
  for (l=1;l<=gv->L;l++) {
    pts[l]=1.0/(2.0*PI*gv->FL[l]);
    eta[l]=gv->DT/pts[l];
  }
  
  ws=2.0*PI/gv->TS;
  
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

  log_infoc(0, "Reading Qp model: %s.qp\n",gv->MFILE);
  sprintf(filename,"%s.qp",gv->MFILE);
  FILE *fp_qp=fopen(filename,"r");
  if (!fp_qp) log_fatal("Could not open Qp model %s.\n", filename);
  
  log_infoc(0, "Reading Qs model: %s.qs\n",gv->MFILE);
  sprintf(filename,"%s.qs",gv->MFILE);
  FILE *fp_qs=fopen(filename,"r");
  if (!fp_qs) log_fatal("Could not open Qs model %s.\n", filename);

  /* loop over global grid */
  for (int i=1;i<=gv->NXG;i++){
    for (int j=1;j<=gv->NYG;j++){
      fread(&vp, sizeof(float), 1, fp_vp);
      fread(&vs, sizeof(float), 1, fp_vs);
      fread(&rhov, sizeof(float), 1, fp_rho);
      fread(&qp, sizeof(float), 1, fp_qp);
      fread(&qs, sizeof(float), 1, fp_qs);

      tp=2.0/(qp*gv->L); 
      ts=2.0/(qs*gv->L);
      muv=vs*vs*rhov;
      piv=vp*vp*rhov;

      sumu=0.0;
      sumpi=0.0;
      for (l=1;l<=gv->L;l++){
	sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
	sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
      }

      muv=muv/(1.0+sumu);
      piv=piv/(1.0+sumpi);
      
      /* only the PE which belongs to the current global gridpoint 
	 is saving model parameters in his local arrays */
      if ((gv->POS[1]==((i-1)/gv->NX)) && 
	  (gv->POS[2]==((j-1)/gv->NY))){
	ii=i-gv->POS[1]*gv->NX;
	jj=j-gv->POS[2]*gv->NY;
	
	taus[jj][ii]=ts;
	taup[jj][ii]=tp;
	u[jj][ii]=muv;
	rho[jj][ii]=rhov;
	pi[jj][ii]=piv;
      }
    }
  }

  fclose(fp_vp);
  fclose(fp_vs);
  fclose(fp_rho);
  fclose(fp_qp);
  fclose(fp_qs);
  
  free_vector(pts,1,gv->L);

  return;
}
