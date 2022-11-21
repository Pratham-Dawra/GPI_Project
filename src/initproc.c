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
/* -------------------------------------------------------------
 * This is function initproc.
 * Dividing the 2-D FD grid into domains and assigning the
 * PEs to these domains,
 *
 * -------------------------------------------------------------*/

#include "fd.h"

void initproc(GlobVar *gv)	{

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	if ((gv->NPROC != gv->NP)  && (MYID==0)) {
		fprintf(gv->FP,"You specified NPROC =  %d (in parameter file) and NP = %d (command line) \n",gv->NPROC,gv->NP);
		declare_error("NP and NPROC differ !");
	}

	/*C-- determine the length of the subarray on this processor*/
	gv->IENDX = gv->NX/gv->NPROCX;
	gv->IENDY = gv->NY/gv->NPROCY;

	/* POS(1) indicates x POSition of the processor in the 
		     logical 3D processor array*/
	if ((gv->NX%gv->NPROCX)>0)
		declare_error(" NX%NPROX (modulus) must be zero  !");
	if ((gv->NY%gv->NPROCY)>0)
		declare_error(" NY%NPROY (modulus) must be zero  !");

	if (MYID==0){
		fprintf(gv->FP,"\n **Message from initprocs (printed by PE %d):\n",MYID);
		fprintf(gv->FP," Size of subarrays in gridpoints:\n");
		fprintf(gv->FP," IENDX= %d\n",gv->IENDX);
		fprintf(gv->FP," IENDY (vertical) = %d\n",gv->IENDY);
	}

	/*---------------   index is indicating neighbouring processes	--------------------*/
	gv->INDEX[1]=MYID-1;           /* left	*/
	gv->INDEX[2]=MYID+1;           /* right	*/
	gv->INDEX[3]=MYID-gv->NPROCX;  /* upper	*/
	gv->INDEX[4]=MYID+gv->NPROCX;  /* lower	*/
	
	/*---------------   POS indicates the processor location in the 3D logical processor array	---------*/
	gv->POS[1] = MYID%gv->NPROCX;  /*  x coordinate */
	gv->POS[2] = (MYID/gv->NPROCX);/*  y coordinate */

	if (gv->POS[1] == 0)        gv->INDEX[1]=gv->INDEX[1] + gv->NPROCX;        	  
	if (gv->POS[1] == gv->NPROCX-1) gv->INDEX[2]=gv->INDEX[2] - gv->NPROCX;          	 
	if (gv->POS[2] == 0)        gv->INDEX[3]=(gv->NPROCX*gv->NPROCY)+MYID-gv->NPROCX; 	 
	if (gv->POS[2] == gv->NPROCY-1) gv->INDEX[4]=MYID+gv->NPROCX-(gv->NPROCX*gv->NPROCY);	 

	fprintf(gv->FP,"\n");
	fprintf(gv->FP," **Message from initprocs (written by PE %d):\n",MYID);
	fprintf(gv->FP," Processor locations in the 2D logical processor array\n");
	fprintf(gv->FP," MYID \t POS(1):left,right \t POS(2): top, bottom\n");

	fprintf(gv->FP," %d \t\t %d: %d,%d \t\t %d: %d,%d \n",
	    MYID,gv->POS[1],gv->INDEX[1],gv->INDEX[2],gv->POS[2],gv->INDEX[3],gv->INDEX[4]);
}
