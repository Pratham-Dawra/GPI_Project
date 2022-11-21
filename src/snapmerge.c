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
 *   loop over snapshotfiles which have to be merged.                                   
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int main(int argc, char **argv) {

  int  nsnap;
  char *fileinp="";
  GlobVar gv = {.FP=stdout};

  /* ============================================== */
  /* Open parameter-file to check if auto mode or not*/
  fileinp = argv[1];

  printf(" ***********************************************************\n");
  printf(" This is program SNAPMERGE.\n");
  printf(" Merge of snapshot files from the parallel\n 2-D Viscoelastic Finite Difference Modeling\n\n");
  printf(" written by  T. Bohlen\n");
  printf(" Geophysical Institute, Department of Physics,\n");
  printf(" Karlsruhe Institute of Technology, Karlsruhe, Germany\n");
  printf(" http://www.gpi.kit.edu\n");
  printf(" ***********************************************************\n\n");
  printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi2D.json\n");
  printf(" Input file for the snapmerge process from command line : %s\n",fileinp);
  
  if ((gv.FP=fopen(fileinp,"r"))==NULL) {
      declare_error(" Opening input file failed.");
  } else {
      printf(" Opening input file was successful.\n\n");
  }
  fclose(gv.FP);

  /* =================================================== */
  /* read standard input file */
  gv.FP = stdout;
  if (strstr(fileinp,".json")) {
      read_par_json(fileinp, &gv);
  } else {
      printf("Parameter file has no .json suffix.");
      exit(EXIT_FAILURE);
  }

  gv.NXG = gv.NX;
  gv.NYG = gv.NY;	
  gv.NX = gv.NXG/gv.NPROCX;
  gv.NY = gv.NYG/gv.NPROCY;
  
  nsnap=1+iround((gv.TSNAP2-gv.TSNAP1)/gv.TSNAPINC);
  
  gv.FP = stdout;
  
  switch(gv.SNAP){
  case 1 : /*particle velocity*/
    merge(nsnap,1,&gv);
    merge(nsnap,2,&gv);
    break;
  case 2 : /*pressure */
    merge(nsnap,6,&gv);
    break;
  case 4 : /*particle velocity*/
    merge(nsnap,1,&gv);
    merge(nsnap,2,&gv);
    merge(nsnap,6,&gv);
  case 3 :
    merge(nsnap,4,&gv);
    merge(nsnap,5,&gv);
    break;
  default :
    warning(" snapmerge: cannot identify content of snapshot !");
    break;
  }
  return EXIT_SUCCESS;	
}
