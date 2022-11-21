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

#include "fd.h"


void mergemod(char modfile[STRING_SIZE], int format, GlobVar *gv){

	char file[STRING_SIZE];
	FILE *fp[gv->NPROCY][gv->NPROCX], *fpout;
	int i, j, ip, jp;
	float a;

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	int remove(const char * filename);

	fprintf(gv->FP,"\n **Message from mergemod (printed by PE %d):\n",MYID);
	fprintf(gv->FP," PE %d starts merge of %d model files \n",MYID,gv->NPROC);	

	fprintf(gv->FP,"\n writing merged model file to  %s \n",modfile);
	fpout=fopen(modfile,"w");

	fprintf(gv->FP," Opening model files: %s.??? ",modfile);
	for (ip=0;ip<=gv->NPROCX-1; ip++)
  	for (jp=0;jp<=gv->NPROCY-1; jp++){
      		sprintf(file,"%s.%i.%i",modfile,ip,jp);
      		fp[jp][ip]=fopen(file,"r");
      		if (fp[jp][ip]==NULL) declare_error("merge: can't read model file !"); 
      	}

	fprintf(gv->FP," ... finished. \n");
    
	fprintf(gv->FP," Copying...");

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
	fprintf(gv->FP," ... finished. \n");

	for (ip=0;ip<=gv->NPROCX-1; ip++)
   	for (jp=0;jp<=gv->NPROCY-1; jp++){
      		fclose(fp[jp][ip]);
    }
	fclose(fpout);
	
	fprintf(gv->FP," Use \n");
	fprintf(gv->FP," ximage n1=%d < %s  label1=Y label2=X title=%s \n",
      			((gv->NYG-1)/gv->IDY)+1,modfile,modfile);
	fprintf(gv->FP," to visualize model. \n");

	fprintf(gv->FP," Removing model files produced by PEs \n");
	for (ip=0;ip<=gv->NPROCX-1; ip++)
  	for (jp=0;jp<=gv->NPROCY-1; jp++){
  		sprintf(file,"%s.%i.%i",modfile,ip,jp);
  		if (remove(file) == 0) {
            fprintf(gv->FP, "The file %s is deleted successfully.",file);
        } else {
            printf("The file %s is not deleted.",file);
        }
	}
}
