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
 *   merge snapshots files written by the different processes to 
 *   a single file                                 
 *
 *  ----------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "fd.h"

void merge(int nsnap, int type, GlobVar *gv){

	char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10];
	FILE *fp[gv->NPROCY][gv->NPROCX], *fpout;
	int i, j, ip, jp, n;
	float a, max;

	switch((gv->SNAP_FORMAT)){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}
    
	switch(type){
	case 1:
        fprintf(gv->FP," x-component of particle velocity");
		strcat(ext,".vx");
		break;
	case 2:
		fprintf(gv->FP," y-component of particle velocity");
		strcat(ext,".vy");
		break;
	case 4:
		fprintf(gv->FP," P-wave energyfield");
		strcat(ext,".div");
		break;
	case 5:
		fprintf(gv->FP," S-wave energyfield");
		strcat(ext,".curl");
		break;
	case 6:
		fprintf(gv->FP," pressure");
		strcat(ext,".p");
		break;
	default:
		declare_error(" merge: cannot find snapfiles! ");
		break;
	}

	sprintf(mfile,"%s%s",(gv->SNAP_FILE),ext);
	fprintf(gv->FP," (files: %s.??? ).\n",mfile);
	sprintf(outfile,"%s%s",(gv->SNAP_FILE),ext);
	fprintf(gv->FP,"\n writing merged snapshot file to  %s\n",outfile);
	fpout=fopen(outfile,"w");
	fprintf(gv->FP," Opening snapshot files: %s.??? ",mfile);

	for (ip=0;ip<=(gv->NPROCX)-1; ip++)
		for (jp=0;jp<=(gv->NPROCY)-1; jp++){
			sprintf(file,"%s.%i.%i",mfile,ip,jp);
			fp[jp][ip]=fopen(file,"r");
			if (fp[jp][ip]==NULL) declare_error("merge: can't read snapfile !");
		}

	fprintf(gv->FP," ... finished.\n Copying...");

	max=0.0;
	for (n=0;n<=nsnap-2; n++){
		for (ip=0;ip<=(gv->NPROCX)-1; ip++){
			for (i=1;i<=(gv->NX);i+=(gv->IDX)){
				for (jp=0;jp<=(gv->NPROCY)-1; jp++){
					for (j=1;j<=(gv->NY);j+=(gv->IDY)){
						a=readdsk(fp[jp][ip],(gv->SNAP_FORMAT));
						if (a>max) max=a;
						writedsk(fpout,a,(gv->SNAP_FORMAT));
					}
				}
			}
		}
	}
	fprintf(gv->FP," ... finished.\n");

	for (ip=0;ip<=(gv->NPROCX)-1; ip++)
		for (jp=0;jp<=(gv->NPROCY)-1; jp++){
			fclose(fp[jp][ip]);
		}

	if ((gv->SNAP_FORMAT)==3){
		fprintf(gv->FP," Use \n");
		fprintf(gv->FP," xmovie n1=%d n2=%d  < %s loop=1 label1=Y label2=X title=%%g d1=%f d2=%f f1=%f f2=%f clip=%e\n",
				(((gv->NYG)-1)/(gv->IDY))+1,(((gv->NXG)-1)/(gv->IDY))+1,outfile,gv->DH,gv->DH,gv->DH,gv->DH,max/10.0);
		fprintf(gv->FP," to play the movie.\n");
		fprintf(gv->FP," ==============================================================\n");
	}
}


