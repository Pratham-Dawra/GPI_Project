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
/*-------------------------------------------------------------
 *  Check FD-Grid for stability and grid dispersion.
 *  If the stability criterion is not fulfilled the program will
 *  terminate.
 *
 *  ----------------------------------------------------------*/

#include <limits.h>

#include "fd.h"
#include "logging.h"

/* TODO: Calculation of global grid size wrong! To be fixed. */

void checkfd(float **prho, float **ppi, float **pu,
             float **ptaus, float **ptaup, float *peta, float *hc, float **srcpos, int nsrc, int **recpos, int ntr, GlobVar *gv) 
{
  float  c, cmax_p=0.0, cmin_p=1e9, cmax_s=0.0, cmin_s=1e9,  cwater=1.0, fmax, gamma;
  float  cmax=0.0, cmin=1e9, sum, dtstab, dhstab, ts, cmax_r, cmin_r, temporal;
  float snapoutx=0.0, snapouty=0.0;
  float srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
  float srec_maxx=-1.0, srec_maxy=-1.0;
  float CFL;
  const float w=2.0*PI/gv->TS; /*center frequency of source*/

  int i, j, k, l, ny1=1, nx, ny, nfw;

  int MYID;
  MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    
  nx=gv->NX;
  ny=gv->NY;

  /* We already perform basic I/O checks for all MPI ranks in check_fs - these additional 
     filesystem checks here seem somewhat superfluous... Deactivated for now. */

/* 	fprintf(gv->FP,"\n **********************************************************"); */
/* 	fprintf(gv->FP,"\n ************ CHECKS OF INPUT FILE PARAMETERS  ************"); */
/* 	fprintf(gv->FP,"\n **********************************************************\n\n"); */
/* 	fprintf(gv->FP,"\n **Message from checkfd_ssg (printed by PE %d):\n",MYID); */
/* 	fprintf(gv->FP,"\n\n ------------------ CHECK OUTPUT FILES --------------------------\n"); */
/* 	/\* The original checks might delete files accidentally that would not be overwritten anyway. */
/* 	   and did not test accessibility from all CPUs which may be vary, especially in distributed clusters *\/ */
/* 	/\*Checking SNAP Output *\/ */
/* 	/\*-------------------- *\/ */
/* 	/\*-------------------------------------- *\/ */
/* 	/\* only MYID is performing the file checks, if to many PEs try to do it simultaneously, */
/* 	 * the code and/or FILESYSTEM and/or MPI implementation will cause segmentation faults */
/* 	 *\/ */
/* 	if ((gv->SNAP>0) && (MYID==0)) { */
/* 		switch (gv->SNAP_FORMAT) { */
/* 			case 1: */
/* 				sprintf(file_ext, ".su"); */
/* 				strcpy(xmod, "ab"); */
/* 				break; */
/* 			case 2: */
/* 				sprintf(file_ext, ".asc"); */
/* 				strcpy(xmod, "a"); */
/* 				break; */
/* 			case 3: */
/* 				sprintf(file_ext, ".bin"); */
/* 				strcpy(xmod, "ab"); */
/* 				break; */
/* 			default: */
/* 				declare_error(" Sorry. Snapshot format (SNAP_FORMAT) unknown. \n"); */
/* 				break; */
/* 		} */
/* 		fprintf(gv->FP," Check accessibility for snapshot files ... \n"); */
/* 		switch (gv->SNAP) { */
/* 			case 1 : */
/* 				sprintf(xfile,"%s%s.vx.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0  cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				sprintf(xfile,"%s%s.vy.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				break; */
/* 			case 2 : */
/* 				sprintf(xfile,"%s%s.p.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				break; */
/* 			case 4 : */
/* 				sprintf(xfile,"%s%s.vx.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				sprintf(xfile,"%s%s.vy.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				sprintf(xfile,"%s%s.p.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 			case 3 : */
/* 				sprintf(xfile,"%s%s.div.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				sprintf(xfile,"%s%s.curl.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]); */
/* 				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile); */
/* 				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile); */
/* 				else fclose(fpcheck); */
/* 				break; */
/* 		} */
/* 	} */
/* 	/\*Checking SEISMOGRAM Output Particle velocities *\/ */
/* 	/\*-------------------------------------- *\/ */
/* 	/\* only MYID is performing the file checks, if to many PEs try to do it simultaneously, */
/* 	 * the code and/or FILESYSTEM and/or MPI implementation will cause segmentation faults */
/* 	 *\/ */
/* 	if ((gv->SEISMO>0) && (MYID==0)) { */
/* 		fprintf(gv->FP," Check accessibility for seismogram files of each PE ... \n"); */
/* 		fprintf(gv->FP," However, the list below refers only to PE0  ... \n"); */
/* 		switch (gv->SEIS_FORMAT) { */
/* 			case 1: */
/* 				sprintf(file_ext,"su"); */
/* 				break; */
/* 			case 2: */
/* 				sprintf(file_ext,"txt"); */
/* 				break; */
/* 			case 3: */
/* 				sprintf(file_ext,"bin"); */
/* 				break; */
/* 		} */
/* 		if (gv->SEIS_FORMAT==2) strcpy(xmod,"a"); */
/* 		else strcpy(xmod,"w+b"); */
/* 		/\*MYID=0 is checking if all seismogram file written by other PEs can be written */
/* 		 * assuming that all PEs can address the files ystem equally */
/* 		 *\/ */
/* 		for (myidcounter=0; myidcounter< (gv->NPROCX*gv->NPROCY); myidcounter++) { */
/* 			switch (gv->SEISMO) { */
/* 				case 1: /\* particle velocities only *\/ */
/* 					sprintf(xfile,"%s_vx.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					/\*in case of number of PE's=500, there will be 500 messages, too many to be displayed! *\/ */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					sprintf(xfile,"%s_vy.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					break; */
/* 				case 2 : /\* pressure only *\/ */
/* 					sprintf(xfile,"%s_p.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					break; */
/* 				case 4 : /\* everything *\/ */
/* 					sprintf(xfile,"%s_vx.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					/\*in case of number of PE's=500, there will be 500 messages, too many to be displayed! *\/ */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					sprintf(xfile,"%s_vy.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					sprintf(xfile,"%s_p.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 				case 3 : /\* curl and div only *\/ */
/* 					sprintf(xfile,"%s_div.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					sprintf(xfile,"%s_curl.%s.%d",gv->SEIS_FILE,file_ext,myidcounter); */
/* 					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile); */
/* 					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile); */
/* 					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage); */
/* 					else fclose(fpcheck); */
/* 					remove(xfile); */
/* 					break; */
/* 			} */
/* 		} */
/* 	} */

/* 	/\*Checking CHECKPOINT Output *\/ */
/* 	/\*-------------------------- *\/ */
/* 	if (gv->CHECKPTREAD>0) { */
/* 		strcpy(xmod,"rb"); */
/* 		sprintf(xfile,"%s.%d",gv->CHECKPTFILE,MYID); */
/* 		fprintf(gv->FP," Check readability for checkpoint files %s... \n",xfile); */
/* 		if (((fpcheck=fopen(xfile,xmod)) ==NULL) && (MYID==0)) declare_error(" PE 0 cannot read checkpoints!"); */
/* 		else fclose(fpcheck); */
/* 	} */
/* 	if ((gv->CHECKPTWRITE>0)) { */
/* 		strcpy(xmod,"ab"); */
/* 		sprintf(xfile,"%s.%d",gv->CHECKPTFILE,MYID); */
/* 		fprintf(gv->FP," Check writability for checkpoint files %s... \n",xfile); */
/* 		if (((fpcheck=fopen(xfile,xmod)) ==NULL) && (MYID==0)) declare_error(" PE 0 cannot write checkpoints!"); */
/* 		else fclose(fpcheck);    /\* Is there any reason to remove it? *\/ */
/* 	} */
/* 	fprintf(gv->FP," Accessibility of output files from PE %d has been checked successfully.\n", MYID); */

   if (0==MYID) {
     log_info("------------------------- Min/max velocities ----------------\n");
     log_info("Note: If any velocity is set <1m/s, it will be ignored while\n");
     log_info("determining stable DH and DT parameters.\n");
   }
  
  /* low Q frame not yet applied as a absorbing boundary */
  /* if (!FREE_SURF) ny1=1+FW;*/
  /*nfw=FW; check only outside the absorbing frame */
  nfw=0;
  cmax_s=0;
  cmin_s=10000;
  cmax_p=0;
  cmin_p=10000;

  /* find maximum model phase velocity of shear waves at infinite frequency within the whole model */
  if (gv->L>0) {   /*viscoelastic*/
    for (i=1+nfw; i<= (nx-nfw); i++) {
      for (j=ny1; j<= (ny-nfw); j++) {
	c=sqrt(pu[j][i]* (1.0/prho[j][i]) * (1.0+gv->L*ptaus[j][i]));
	/*if c is close to zero (water, air), c will be ignored for finding cmax,cmin*/
	if ((cmax_s<c) && (c>cwater)) cmax_s=c;
	/* find minimum model phase velocity of shear waves at center frequency of the source */
	sum=0.0;
	for (l=1; l<=gv->L; l++) {
	  ts=gv->DT/peta[l];
	  sum=sum+ ((w*w*ptaus[j][i]*ts*ts) / (1.0+w*w*ts*ts));
	}
	c=sqrt((pu[j][i]/prho[j][i]));
	if ((cmin_s>c) && (c>cwater)) cmin_s=c;
      }
    }
  } else { /* L=0, elastic */
    for (i=1+nfw; i<= (nx-nfw); i++) {
      for (j=ny1; j<= (ny-nfw); j++) {
	c=sqrt(pu[j][i]/prho[j][i]);
	
	/*if c is close to zero (water, air), c will be ignored for finding	cmax,cmin*/
	if ((c>cwater) && (cmax_s<c)) cmax_s=c;
	if ((c>cwater) && (cmin_s>c)) cmin_s=c;
      }
    }
  }

  /* find maximum model phase velocity of P-waves at infinite frequency within the whole model */
  if (gv->L>0) {   /*viscoelastic*/
    for (i=1+nfw; i<= (nx-nfw); i++) {
      for (j=ny1; j<= (ny-nfw); j++) {
	c=sqrt(ppi[j][i]* (1.0/prho[j][i]) * (1.0+gv->L*ptaup[j][i]));
	if ((c>cwater) && (cmax_p<c)) cmax_p=c;
	/* find minimum model phase velocity of P-waves at center frequency of the source */
	sum=0.0;
	for (l=1; l<=gv->L; l++) {
	  ts=gv->DT/peta[l];
	  sum=sum+ ((w*w*ptaup[j][i]*ts*ts) / (1.0+w*w*ts*ts));
	}
	c=sqrt((ppi[j][i]/prho[j][i]));
	if ((c>cwater) && (cmin_p>c)) cmin_p=c;
      }
    }
  } else { /* L=0, elastic */
    for (i=1+nfw; i<= (nx-nfw); i++) {
      for (j=ny1; j<= (ny-nfw); j++) {
	c=sqrt(ppi[j][i]/prho[j][i]);
	/*if c is close to zero (water, air), c will be ignored for finding	cmax,cmin*/
	if ((c>cwater) && (cmax_p<c)) cmax_p=c;
	if ((c>cwater) && (cmin_p>c)) cmin_p=c;
      }
    }
  }

  log_debug("Vp_min(f=fc)=%.2f, Vp_max(f=inf)=%.2f, Vs_min(f=fc)=%.2f, Vs_max(f=inf)=%.2f\n", cmin_p, cmax_p, cmin_s, cmax_s);
	    
  if (cmax_s>cmax_p) cmax=cmax_s;
  else cmax=cmax_p;

  if (cmin_s<cmin_p) cmin=cmin_s;
  else cmin=cmin_p;

  /* find global maximum for Vp and global minimum for Vs*/
  MPI_Allreduce(&cmax,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&cmin,&cmin_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
  cmax=cmax_r;
  cmin=cmin_r;

  if (gv->FDORDER_TIME==4) {
    temporal=3.0/2.0;
  } else {
    temporal=1.0;
  }

  fmax=2.0/gv->TS;
  dhstab = (cmin/ (hc[0]*fmax));
  gamma = fabs(hc[1]) + fabs(hc[2]) + fabs(hc[3]) + fabs(hc[4]) + fabs(hc[5]) + fabs(hc[6]);
  dtstab = gv->DH/ (sqrt(2) *gamma*cmax*temporal);
  CFL=cmax*gv->DT/gv->DH;

  if (MYID == 0) {
    log_info("Min and max model velocity: V_min=%.2fm/s, V_max=%.2fm/s\n", cmin, cmax);
    log_info("------------------------- Grid dispersion -------------------\n");
    log_info("To limit grid dispersion, the number of grid points per min.\n");
    log_info("wavelength (of S-waves) should be at least 6. Here, the min.\n");
    log_info("wavelength is assumed to be the minimum model phase velocity\n");
    log_info("(of S-waves) at max. frequency divided by the max. frequency\n");
    log_info("of the source (here approximately %.3fHz). The minimum\n",  2.0/gv->TS); 
    log_info("wavelength in the following simulation will be %.3fm. Thus,\n", cmin/fmax);
    log_info("the recommended value is DH=%.3fm. Your value: DH=%.3fm\n", dhstab, gv->DH);
    if (gv->DH>dhstab) {
      log_warn("Grid dispersion will influence wave propagation, choose a smaller grid spacing (DH).");
    }

    log_info("------------------------- Stability check -------------------\n");
    log_info("The simulation is stable if p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n");
    log_info("where cmax is the maximum phase velocity at infinite frequency\n");
    log_info("and gamma=sum(|FD coeff.|). In this simulation, cmax=%.2fm/s.\n", cmax);
    log_info("DT is the time step and DH is the grid size. The Courant-\n");
    log_info("Friedrichs-Lewy number is %.4f. The stability limit for the\n", CFL);
    log_info("time step is DT=%es. Your value: DT=%es\n", dtstab, gv->DT);
    if (gv->DT>dtstab) {
      log_error("The simulation will get unstable, choose a smaller DT.\n");
      log_fatal("Instability of the simulation, update your parameters.\n");
    } else {
      log_info("The simulation will be stable.\n");
    }

    log_info("------------------------- Additional checks -----------------\n");
    if (gv->SNAP) {
      if (gv->TSNAP2>gv->TIME) {
	log_warn("TSNAP2=%e (last snapshot) > time of wave propagation %e. TSNAP2 changed to be equal to TIME.\n",gv->TSNAP2, gv->TIME);
	gv->TSNAP2=gv->TIME;
      }
      snapoutx=gv->NX/ (float) gv->IDX;
      snapouty=gv->NY/ (float) gv->IDY;
      log_info("Output of snapshot grid points (x) per node: %8.2f\n", snapoutx);
      log_info("Output of snapshot grid points (y) per node: %8.2f\n", snapouty);
      if (snapoutx- (int) snapoutx>0) log_fatal("Ratio NX-NPROCX-IDX must be integer.\n");
      if (snapouty- (int) snapouty>0) log_fatal("Ratio NY-NPROCY-IDY must be integer.\n");
    } else {
      log_info("Skipping checks of snapshot parameters.\n");
    }
    
    if (gv->SEISMO) {
      log_info("Number of modeling time steps: %d\n", gv->NT);
      log_info("Seismogram sampling interval in time steps: %d\n", gv->NDT);
      log_info("Number of seismogram output samples: %d\n", gv->NT/gv->NDT);
      
      /* SU and SEG-Y allow 32767 samples, furthermore the exist programs allow for 65535 
	 samples and pseudo-SEG-Y formats allow foralmost arbitrarily long traces.
	 For binary and textual output the limit is arbitrary. USHRT_MAX is the limut of 
	 an unsigned short specified in limits.h */
      
      if ((gv->SEIS_FORMAT==1) && (gv->NT/gv->NDT) > (USHRT_MAX)) {
	log_error("Maximum number of samples per trace in SU format: %d. Your value: %d\n", USHRT_MAX, gv->NT/gv->NDT);
	log_fatal("Too many output samples per receiver for SU format.\n");
      }

      srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
      srec_maxx=-1.0, srec_maxy=-1.0;
      log_info("Checking for receiver position(s) as specified in json file.\n");
      log_info("Global grid size in m: %5.2f (x) : %5.2f (y)\n",gv->NX*gv->DH*gv->NPROCX,gv->NY*gv->DH*gv->NPROCY);    
      if (gv->FREE_SURF==0) log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n", 
				     (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH, (float) gv->FW*gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);     
      if (gv->FREE_SURF==1) log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n", 
				     (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH,gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);
      
      /* find maximum and minimum source positions coordinate ---- from input file*/    
      if (gv->READREC==0) {
	if (gv->XREC1>gv->XREC2) {
	  srec_maxx=gv->XREC1;
	  srec_minx=gv->XREC2;
	} else {
	  srec_maxx=gv->XREC2;
	  srec_minx=gv->XREC1;
	}
	if (gv->YREC1>gv->YREC2) {
	  srec_maxy=gv->YREC1;
	  srec_miny=gv->YREC2;
	} else {
	  srec_maxy=gv->YREC2;
	  srec_miny=gv->YREC1;
	}
	log_info("Number of receiver positions: %d\n", ntr);
      }
      
      if (gv->READREC==1) {
	/* find maximum and minimum source positions coordinate ---- from receiver file*/
	for (k=1; k<=ntr; k++) {
	  /* find maximum source positions coordinate*/
	  if ((recpos[1][k]*gv->DH) >srec_maxx) srec_maxx=recpos[1][k]*gv->DH;
	  if ((recpos[2][k]*gv->DH) >srec_maxy) srec_maxy=recpos[2][k]*gv->DH;
	  /* find minimum source positions coordinate*/
	  if ((recpos[1][k]*gv->DH) <srec_minx) srec_minx=recpos[1][k]*gv->DH;
	  if ((recpos[2][k]*gv->DH) <srec_miny) srec_miny=recpos[2][k]*gv->DH;
	}
	log_info("Number of receiver positions: %d\n", ntr);
      }
      
      log_info("Minimum receiver position coordinates: %5.2f (x) : %5.2f (y)\n",srec_minx,srec_miny);
      log_info("Maximum receiver position coordinates: %5.2f (x) : %5.2f (y)\n",srec_maxx,srec_maxy);

      /* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
      if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
	log_fatal("Coordinate of at least one receiver location is outside the global grid.\n");
      }
      if ((srec_maxx>gv->NX*gv->DH*gv->NPROCX) || (srec_maxy>gv->NY*gv->DH*gv->NPROCY)) {
	log_fatal("Coordinate of at least one receiver location is outside the global grid.\n");
      }
      
      /* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
      if ((srec_maxx< ((float) gv->FW*gv->DH)) || (srec_minx< ((float) gv->FW*gv->DH))) {
	/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
	log_warn("Coordinate of at least one receiver location is inside the Absorbing Boundary (left boundary).\n");
      }
      
      if (srec_maxx> (gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one receiver location is inside the Absorbing Boundary (right boundary).\n");
      }
      
      if (srec_maxy> (gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one receiver location is inside the Absorbing Boundary (lower boundary).\n");
      }
      
      if ((srec_miny< ((float) gv->FW*gv->DH)) && !(gv->FREE_SURF)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one receiver location is inside the Absorbing Boundary (top boundary).\n");
      }
    } else {
      log_info("Skipping checks of seismogram parameters.\n");
    }
   
    if (gv->SRCREC==1) {
      srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
      srec_maxx=-1.0, srec_maxy=-1.0;
      log_info("Checking for source position(s) as specified in source file.\n");
      log_info("Global grid size in m: %5.2f (x) : %5.2f (y)\n",gv->NX*gv->DH*gv->NPROCX,gv->NY*gv->DH*gv->NPROCY);
      if (gv->FREE_SURF==0) log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n", 
				     (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH, (float) gv->FW*gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);
      if (gv->FREE_SURF==1) log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n", 
				     (float) gv->FW*gv->DH,gv->NX* (float) gv->DH*gv->NPROCX- (float) gv->FW*gv->DH,gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);
      
      for (k=1; k<=nsrc; k++) {
	/* find maximum source positions coordinate*/
	if (srcpos[1][k]>srec_maxx) srec_maxx=srcpos[1][k];
	if (srcpos[2][k]>srec_maxy) srec_maxy=srcpos[2][k];
	/* find minimum source positions coordinate*/
	if (srcpos[1][k]<srec_minx) srec_minx=srcpos[1][k];
	if (srcpos[2][k]<srec_miny) srec_miny=srcpos[2][k];
      }

      log_info("Number of source positions: %d\n", nsrc);
      log_info("Minimum source position coordinates: %5.2f (x) : %5.2f (y)\n",srec_minx,srec_miny);
      log_info("Maximum source position coordinates: %5.2f (x) : %5.2f (y)\n",srec_maxx,srec_maxy);

      /* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
      if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
	log_fatal("Coordinate of at least one source location is outside the global grid.\n");
      }
      
      if ((srec_maxx>gv->NX*gv->DH*gv->NPROCX) || (srec_maxy>gv->NY*gv->DH*gv->NPROCY)) {
	log_fatal("Coordinate of at least one source location is outside the global grid.\n");
      }
      
      /* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
      if ((srec_maxx< ((float) gv->FW*gv->DH)) || (srec_minx< ((float) gv->FW*gv->DH))) {
	/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
	log_warn("Coordinate of at least one source location is inside the Absorbing Boundary (left boundary).\n");
      }
      if (srec_maxx> (gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one source location is inside the Absorbing Boundary (right boundary).\n");
      }
      if (srec_maxy> (gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one source location is inside the Absorbing Boundary (lower boundary).\n");
      }
      if ((srec_miny< ((float) gv->FW*gv->DH)) && !(gv->FREE_SURF)) {
	/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
	log_warn("Coordinate of at least one source location is inside the Absorbing Boundary (top boundary).\n");
      }
    }

    log_info("------------------------- Absorbing frame checks ------------\n");
    if (gv->ABS_TYPE==1) {
      log_info("CPML boundary (ABS_TYPE=1) with width of %d grid points (%5.2fm).\n",gv->FW, (float) gv->FW*gv->DH);
      if (gv->FW<10) {
	log_warn("Width (FW) of absorbing frame should be at least 10 grid points.\n");
	log_warn("Be aware of artificial reflections from grid boundaries!\n");
      }
    }
    if (gv->ABS_TYPE==2) {
      log_info("Absorbing boundary (ABS_TYPE=2) with width of %d grid points (%5.2fm).\n",gv->FW, (float) gv->FW*gv->DH);
      if (gv->FW<30) {
	log_warn("Width (FW) of absorbing frame should be at least 30 grid points.\n");
	log_warn("Be aware of artificial reflections from grid boundaries!\n");
      }
    }
    if (((gv->NX) < gv->FW) || ((gv->NY) <gv->FW))	{
      log_error("Width of absorbing boundary (FW=%d grid points) is larger than at least one subdomain dimension.\n",gv->FW);
      log_error("Subdomain dimensions: NX/NPROCX=%d, NY/NPROCY=%d in grid points.\n",gv->NX,gv->NY);
      log_fatal("You need to choose a smaller width of absorbing frame (FW) or increase the subdomain dimensions.\n");
    }
    if (gv->BOUNDARY) {
      if (gv->ABS_TYPE==1 || gv->ABS_TYPE==2) {
	log_warn("You have activated a periodic boundary and set an absorbing boundary at the same time!\n");
      }
    }
  }
}
