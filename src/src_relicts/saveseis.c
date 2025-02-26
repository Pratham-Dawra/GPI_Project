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
 *   write seismograms of each PE individually to files
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
		float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc,
		int ntr, float ** srcpos_loc, int ishot,int ns, GlobVar *gv) {

	char vxf[STRING_SIZE*2], vyf[STRING_SIZE*2], curlf[STRING_SIZE*2], divf[STRING_SIZE*2], pf[STRING_SIZE*2],file_ext[5];
	int nsrc=1;
    
    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	if (MYID==0) fprintf(fp,"\n **Message from function saveseis (written by PE %d): \n",MYID);

	switch (gv->SEIS_FORMAT){
	case 1: sprintf(file_ext,"su");  break;
	case 2: sprintf(file_ext,"txt"); break;
	case 3: sprintf(file_ext,"bin"); break;
	}

	if (gv->RUN_MULTIPLE_SHOTS){

		sprintf(vxf,"%s_vx.%s.shot%d.%d",gv->SEIS_FILE,file_ext,ishot,MYID);
		sprintf(vyf,"%s_vy.%s.shot%d.%d",gv->SEIS_FILE,file_ext,ishot,MYID);
		sprintf(curlf,"%s_curl.%s.shot%d.%d",gv->SEIS_FILE,file_ext,ishot,MYID);
		sprintf(divf,"%s_div.%s.shot%d.%d",gv->SEIS_FILE,file_ext,ishot,MYID);
		sprintf(pf,"%s_p.%s.shot%d.%d",gv->SEIS_FILE,file_ext,ishot,MYID);
	}
	else{

		sprintf(vxf,"%s_vx.%s.%d",gv->SEIS_FILE,file_ext,MYID);
		sprintf(vyf,"%s_vy.%s.%d",gv->SEIS_FILE,file_ext,MYID);
		sprintf(curlf,"%s_curl.%s.%d",gv->SEIS_FILE,file_ext,MYID);
		sprintf(divf,"%s_div.%s.%d",gv->SEIS_FILE,file_ext,MYID);
		sprintf(pf,"%s_p.%s.%d",gv->SEIS_FILE,file_ext,MYID);

	}


	switch (gv->SEISMO){
	case 1 : /* particle velocities only */

		fprintf(fp,"\n PE %d is writing %d seismogram traces (vx)   to %s ",MYID,ntr,vxf);
		outseis(fp,fopen(vxf,"w"),sectionvx,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);
		fprintf(fp,"\n PE %d is writing %d seismogram traces (vy)   to %s ",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),sectionvy,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);

		break;
	case 2 : /* pressure only */
		fprintf(fp,"\n PE %d is writing %d seismogram traces (p)    to %s \n",MYID,ntr,pf);
		outseis(fp,fopen(pf,"w"),sectionp,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);

		break;
	case 3 : /* curl and div only */

		fprintf(fp,"\n PE %d is writing %d seismogram traces (div)  to %s ",MYID,ntr,divf);
		outseis(fp,fopen(divf,"w"),sectiondiv,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);
		fprintf(fp,"\n PE %d is writing %d seismogram traces (curl) to %s ",MYID,ntr,curlf);
		outseis(fp,fopen(curlf,"w"),sectioncurl,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);

		break;
	case 4 : /* everything */
		fprintf(fp,"\n PE %d is writing %d seismogram traces (vx)   to %s ",MYID,ntr,vxf);
		outseis(fp,fopen(vxf,"w"),sectionvx,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);
		fprintf(fp,"\n PE %d is writing %d seismogram traces (vy)   to %s ",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),sectionvy,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);

		fprintf(fp,"\n PE %d is writing %d seismogram traces (p)    to %s ",MYID,ntr,pf);
		outseis(fp,fopen(pf,"w"),sectionp,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);

		fprintf(fp,"\n PE %d is writing %d seismogram traces (div)  to %s ",MYID,ntr,divf);
		outseis(fp,fopen(divf,"w"),sectiondiv,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);
		fprintf(fp,"\n PE %d is writing %d seismogram traces (curl) to %s \n",MYID,ntr,curlf);
		outseis(fp,fopen(curlf,"w"),sectioncurl,recpos,recpos_loc,ntr,srcpos_loc,nsrc,ns,gv->SEIS_FORMAT, ishot, gv);
		break;

	}
}
