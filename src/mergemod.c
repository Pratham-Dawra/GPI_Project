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
 *   merge model files written by the different processes to 
 *   a single file                                 
 *
 *  ----------------------------------------------------------------------*/

#include <stdio.h>
#include "fd.h"
#include "logging.h"

void mergemod(const char* modfile, int format, GlobVar *gv)
{
  char file[STRING_SIZE];
  FILE *fp[gv->NPROCY][gv->NPROCX], *fpout = NULL;
  int i, j, ip, jp;
  float a;

  log_info("Merging model files. Output file: %s\n",modfile);

  fpout=fopen(modfile,"w");
  if (!fpout) log_fatal("Could not open output file %s for writing.\n", modfile);

  log_debug("Opening model files %s.*.* for reading.\n", modfile);

  for (ip=0;ip<=gv->NPROCX-1; ip++) {
    for (jp=0;jp<=gv->NPROCY-1; jp++) {
      sprintf(file,"%s.%i.%i",modfile,ip,jp);
      fp[jp][ip]=fopen(file,"r");
      if (!fp[jp][ip]) log_fatal("Could not open file %s for reading.\n", file); 
    }
  }

  log_debug("Finished opening all model files. Now copying the data.\n");
    
  for (ip=0;ip<=gv->NPROCX-1; ip++){
    for (i=1;i<=gv->NX;i+=gv->IDX){
      for (jp=0;jp<=gv->NPROCY-1; jp++){
	for (j=1;j<=gv->NY;j+=gv->IDY){
	  a=readdsk(fp[jp][ip],format);
	  writedsk(fpout,a,format);
	}
      }
    }
  }
  log_debug("Finished copying the data. Now closing all input files.\n");

  for (ip=0;ip<=gv->NPROCX-1; ip++) {
    for (jp=0;jp<=gv->NPROCY-1; jp++) {
      fclose(fp[jp][ip]);
    }
  }
  log_debug("Finished closing all input files.\n");

  fclose(fpout);

  log_info("Use...\n");
  log_info("%sximage n1=%d < %s label1=\"Y\" label2=\"X\" title=\"%s\"%s\n",
	   LOG_COLOR_BOLD,((gv->NYG-1)/gv->IDY)+1,modfile,modfile,LOG_ALL_RESET);
  log_info("...to visualize the model.\n");
  
  log_debug("Removing individual model files procudes by PEs.\n");
  for (ip=0;ip<=gv->NPROCX-1; ip++) {
    for (jp=0;jp<=gv->NPROCY-1; jp++) {
      sprintf(file,"%s.%i.%i",modfile,ip,jp);
      if (remove(file)!=0) { 
	log_warn("Could not remove model file %s. Continuing...\n",file);
      }
    }
  }
}
