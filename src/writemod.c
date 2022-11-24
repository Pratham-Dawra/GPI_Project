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
 *   write local model to file              
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void writemod(char modfile[STRING_SIZE], float ** array, int format, GlobVar *gv)
{
  int i, j;
  FILE *fpmod = NULL;
  char file[STRING_SIZE];

  sprintf(file,"%s.%i.%i",modfile,gv->POS[1],gv->POS[2]);

  log_debug("Writing model to file %s.\n", file);
	
  fpmod=fopen(file,"w");
  if (!fpmod) log_fatal("Could not open model output file %s for writing.\n", file);
  for (i=1;i<=gv->NX;i+=gv->IDX)
    for (j=1;j<=gv->NY;j+=gv->IDY)
      writedsk(fpmod,array[j][i],format);

  fclose(fpmod);
}


