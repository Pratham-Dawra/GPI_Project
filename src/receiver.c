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
 *  compute receiver positions or read them from external file
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include <stdbool.h>

int **receiver(int *ntr, GlobVar *gv) 
{
  int **recpos1=NULL, **recpos=NULL;
  int   itr=1, itr1=0, itr2=0, recflag=0, n, i, j;
  float nxrec=0, nyrec=0;
  float nxrec1, nxrec2, nyrec1, nyrec2;
  float xrec, yrec;
  bool testbuff1, testbuff2, testbuff3;
  char bufferstring[10], buffer[STRING_SIZE];
  FILE *fpr = NULL;
  
  if (gv->MPID==0) {
    log_info("------------------------- Receiver positions ----------------\n");
    if (gv->READREC) { /* read receiver positions from file */
      log_info("Reading receiver positions from file %s.\n",gv->REC_FILE);
      fpr=fopen(gv->REC_FILE,"r");
      if (fpr==NULL) log_fatal("Receiver file %s could not be opened.\n",gv->REC_FILE);
      *ntr=0;
      
      /* counts the number of receivers in the receiver file */
      while(fgets(buffer, STRING_SIZE, fpr)){
	testbuff1=strchr(buffer,'#');
	testbuff2=strchr(buffer,'%');
	testbuff3=sscanf(buffer,"%s",bufferstring)==1;
	
	/*the following output is for debugging*/
	/*testbuff4=(testbuff1==1 || testbuff2==1);
	  log_std(fp," buffer: _%s_with testbuff1=_%i_ testbuff2=_%i_testbuff3=_%i_ testbuff4=_%i_\n",buffer,testbuff1, testbuff2, testbuff3,testbuff4);*/

	/* checks if the line contains a '%' or '#' character which indicates a
	   comment line, and if the reading of a string was successful, 
	   which is not the case for an empty line*/
	if (((testbuff1==1 || testbuff2==1)==0) && testbuff3==1) ++(*ntr);
      }
      
      rewind(fpr);
      
      recpos1=imatrix(1,3,1,*ntr);
      for (itr=1;itr<=*ntr;itr++){
	fscanf(fpr,"%f%f\n",&xrec, &yrec);
	recpos1[1][itr]=iround((xrec+gv->REFREC[1])/gv->DH);
	recpos1[2][itr]=iround((yrec+gv->REFREC[2])/gv->DH);
	recpos1[3][itr]=itr;
      }
      fclose(fpr);

      /* check if more than one receiver is located at the same gridpoint */
      for (itr=1;itr<=(*ntr-1);itr++)
	for (itr1=itr+1;itr1<=*ntr;itr1++)
	  if ((recpos1[1][itr]==recpos1[1][itr1])
	      && (recpos1[2][itr]==recpos1[2][itr1]))
	    recpos1[1][itr1]=-(++recflag);
      
      recpos=imatrix(1,3,1,*ntr-recflag);
      for (itr=1;itr<=*ntr;itr++)
	if (recpos1[1][itr]>0){
	  recpos[1][++itr2]=recpos1[1][itr];
	  recpos[2][itr2]=recpos1[2][itr];
	  recpos[3][itr2]=recpos1[3][itr];
	}
      
      int ntr_orig = *ntr;
      *ntr=itr2;
      if (recflag>0){
	log_warn("Co-located receivers positions; reduced no. of receivers from %d to %d.\n", ntr_orig, *ntr);
      }
      
      free_imatrix(recpos1,1,3,1,*ntr);  
    }
    else if (gv->REC_ARRAY>0){
      log_info("Generating receiver planes as specified in input file.\n");
      
      *ntr=(1+(gv->NXG-2*gv->FW)/gv->DRX)*gv->REC_ARRAY;
      recpos=imatrix(1,3,1,*ntr);
      itr=0;
      for (n=0;n<=gv->REC_ARRAY-1;n++){
	j=iround((gv->REC_ARRAY_DEPTH+gv->REC_ARRAY_DIST*(float)n)/gv->DH);
	for (i=gv->FW;i<=gv->NXG-gv->FW;i+=gv->DRX){
	  itr++;
	  recpos[1][itr]=i;
	  recpos[2][itr]=j;
	  recpos[3][itr]=itr;
	}
      }
    }
    else {         /* straight horizontal or vertical line of receivers */
      log_info("Receiver positions specified by json parameter file.\n");
      nxrec1=gv->XREC1/gv->DH;   /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
      nyrec1=gv->YREC1/gv->DH;   /* the positions of the first and last receiver*/
      nxrec2=gv->XREC2/gv->DH;	 /* in gridpoints */
      nyrec2=gv->YREC2/gv->DH;
      
      /* only 1 receiver */
      if (((iround(nxrec2)-iround(nxrec1))==0) && ((iround(nyrec2)-iround(nyrec1))==0)) {
	log_info("A single receiver position determined from json file.\n");
	*ntr = 1;
	recpos = imatrix(1,3,1,*ntr);
	recpos[1][1] = iround(nxrec1);
	recpos[2][1] = iround(nyrec1);
	recpos[3][1] = 1;
      } else if (((iround(nyrec2)-iround(nyrec1))==0)) {
	log_info("A horizontal receiver line (in x-direction) determined from json file.\n");
	*ntr=iround((nxrec2-nxrec1)/gv->NGEOPH)+1;
	recpos=imatrix(1,3,1,*ntr);
	nxrec = nxrec1;
	for (n=1;n<=*ntr;n++){
	  nyrec=nyrec1+((nyrec2-nyrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
	  itr=iround((nxrec-nxrec1)/gv->NGEOPH)+1;
	  recpos[1][itr] = iround(nxrec);
	  recpos[2][itr] = iround(nyrec);
	  recpos[3][itr] = itr;
	  nxrec += gv->NGEOPH;
	}
      }
      else if (((iround(nxrec2)-iround(nxrec1))==0)) {
	log_info("A vertical receiver line (in y-direction) determined from json file.\n");
	*ntr=iround((nyrec2-nyrec1)/gv->NGEOPH)+1;
	recpos=imatrix(1,3,1,*ntr);
	nyrec = nyrec1;
	for (n=1;n<=*ntr;n++){
	  nxrec=nxrec1+((nxrec2-nxrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
	  itr=iround((nyrec-nyrec1)/gv->NGEOPH)+1;
	  recpos[1][itr] = iround(nxrec);
	  recpos[2][itr] = iround(nyrec);
	  recpos[3][itr] = itr;
	  nyrec += gv->NGEOPH;
	}
      }
      else {
	/* arbitrary geophone-line */
	log_error("No horizontal or vertical receiver line is specified in the json file.\n");
	log_error("In order to define an arbitrary receiver line, please make use of an external receiver file (READREC=1).\n");
	log_fatal("Error in specifying receiver coordinates in the json file!\n");
      }
    } /* end of if receivers specified in input file */
    log_info("Number of receiver positions: %d\n",*ntr); 
    
    for (itr=1;itr<=*ntr;itr++) {
      /* check for 0's */
      for (n=1;n<=2;n++) {
	if (recpos[n][itr]==0) recpos[n][itr] = 1;
      }
    }
  } /* End of if(gv->MPID==0) */
  

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);
  if (gv->MPID!=0) recpos=imatrix(1,3,1,*ntr);
  MPI_Bcast(&recpos[1][1],(*ntr)*3,MPI_INT,0,MPI_COMM_WORLD);

  if (gv->MPID==0) {
    if (*ntr>50) log_warn("The following table is quite large (%d lines); only printing the first 50 entries!\n",*ntr);
    log_info("Receiver positions in the global model-system:\n");
    log_info(" x (gridpoints) y (gridpoints) \t x (in m) \t y (in m)\n");
    log_info(" -------------  -------------- \t ---------\t --------\n");
    
    int maxprint = *ntr;
    if (maxprint>50) maxprint = 50;

    for (itr=1;itr<=maxprint;itr++) {
      log_info(" %d\t\t %d \t\t %6.2f \t %6.2f\n",recpos[1][itr],recpos[2][itr],recpos[1][itr]*gv->DH,recpos[2][itr]*gv->DH);
    }
  }
  
  return recpos;
}
