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
#include "logging.h"

void readmod_elastic(float **rho, float **pi, float **u, GlobVar *gv) 
{
  
  /* local variables */
  float rhov, muv, piv, vp, vs;
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

  /* loop over global grid */
  for (int i=1;i<=gv->NXG;i++){
    for (int j=1;j<=gv->NYG;j++){
      fread(&vp, sizeof(float), 1, fp_vp);
      fread(&vs, sizeof(float), 1, fp_vs);
      fread(&rhov, sizeof(float), 1, fp_rho);
      muv=vs*vs*rhov;
      piv=vp*vp*rhov;
      /* only the PE which belongs to the current global gridpoint 
	 is saving model parameters in his local arrays */
      if ((gv->POS[1]==((i-1)/gv->NX)) && (gv->POS[2]==((j-1)/gv->NY))) {
	ii=i-gv->POS[1]*gv->NX;
	jj=j-gv->POS[2]*gv->NY;
	u[jj][ii]=muv;
	rho[jj][ii]=rhov;
	pi[jj][ii]=piv;
      }
    }
  }

  fclose(fp_vp);
  fclose(fp_vs);
  fclose(fp_rho);

  return;
}
