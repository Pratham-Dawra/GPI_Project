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


void writemod(char modfile[STRING_SIZE], float ** array, int format, GlobVar *gv){

	int i, j;
	FILE *fpmod;
	char file[STRING_SIZE];

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	fprintf(gv->FP,"\n\n PE %d is writing model to \n",MYID);
	sprintf(file,"%s.%i.%i",modfile,gv->POS[1],gv->POS[2]);
	fprintf(gv->FP,"\t%s\n\n", file);
	fpmod=fopen(file,"w");
	for (i=1;i<=gv->NX;i+=gv->IDX)
	for (j=1;j<=gv->NY;j+=gv->IDY)
		writedsk(fpmod,array[j][i],format);

	fclose(fpmod);
}


