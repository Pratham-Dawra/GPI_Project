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

void checkfd(float **prho, float **ppi, float **pu,
             float **ptaus, float **ptaup, float *peta, float *hc, float **srcpos, int nsrc, int **recpos, int ntr, GlobVar *gv) {


	//extern float DH, DT, TS, TIME, TSNAP2;
	//extern float XREC1, XREC2, YREC1, YREC2;
	//extern int NX, NY, L, MYID, IDX, IDY, NT, NDT, RSG;
	//extern int READREC, NPROCX,NPROCY, SRCREC, FREE_SURF, ABS_TYPE, FW, BOUNDARY,FDORDER_TIME;
	//extern int SNAP, SEISMO, CHECKPTREAD, CHECKPTWRITE, SEIS_FORMAT[6], SNAP_FORMAT, POS[4];
	//extern char SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE], SNAP_FILE[STRING_SIZE];
	//extern char SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];

	/* local variables */

	float  c, cmax_p=0.0, cmin_p=1e9, cmax_s=0.0, cmin_s=1e9,  cwater=1.0, fmax, gamma;
	float  cmax=0.0, cmin=1e9, sum, dtstab, dhstab, ts, cmax_r, cmin_r, temporal;
	float snapoutx=0.0, snapouty=0.0;
	float srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
	float srec_maxx=-1.0, srec_maxy=-1.0;
	float CFL;
	const float w=2.0*PI/gv->TS; /*center frequency of source*/

	int i, j, k, l, ny1=1, nx, ny, myidcounter, nfw;

	char xfile[STRING_SIZE*2], errormessage[STRING_SIZE*2], xmod[4], file_ext[8];
	FILE *fpcheck;

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    
	nx=gv->NX;
	ny=gv->NY;

	fprintf(gv->FP,"\n **********************************************************");
	fprintf(gv->FP,"\n ************ CHECKS OF INPUT FILE PARAMETERS  ************");
	fprintf(gv->FP,"\n **********************************************************\n\n");
	fprintf(gv->FP,"\n **Message from checkfd_ssg (printed by PE %d):\n",MYID);
	fprintf(gv->FP,"\n\n ------------------ CHECK OUTPUT FILES --------------------------\n");
	/* The original checks might delete files accidentally that would not be overwritten anyway.
	   and did not test accessibility from all CPUs which may be vary, especially in distributed clusters */

	/*Checking SNAP Output */
	/*-------------------- */
	/*-------------------------------------- */
	/* only MYID is performing the file checks, if to many PEs try to do it simultaneously,
	 * the code and/or FILESYSTEM and/or MPI implementation will cause segmentation faults
	 */
	if ((gv->SNAP>0) && (MYID==0)) {

		switch (gv->SNAP_FORMAT) {
			case 1:
				sprintf(file_ext, ".su");
				strcpy(xmod, "ab");
				break;

			case 2:
				sprintf(file_ext, ".asc");
				strcpy(xmod, "a");
				break;

			case 3:
				sprintf(file_ext, ".bin");
				strcpy(xmod, "ab");
				break;

			default:
				declare_error(" Sorry. Snapshot format (SNAP_FORMAT) unknown. \n");
				break;
		}


		fprintf(gv->FP," Check accessibility for snapshot files ... \n");

		switch (gv->SNAP) {
			case 1 :
				sprintf(xfile,"%s%s.vx.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0  cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				sprintf(xfile,"%s%s.vy.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				break;

			case 2 :
				sprintf(xfile,"%s%s.p.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				break;

			case 4 :
				sprintf(xfile,"%s%s.vx.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				sprintf(xfile,"%s%s.vy.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				sprintf(xfile,"%s%s.p.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

			case 3 :
				sprintf(xfile,"%s%s.div.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				sprintf(xfile,"%s%s.curl.%i.%i",gv->SNAP_FILE,file_ext,gv->POS[1],gv->POS[2]);
				fprintf(gv->FP,"    Check accessibility for snapshot file %s... \n",xfile);

				if ((fpcheck=fopen(xfile,xmod)) ==NULL) err2(" PE0 cannot write snapshots to %s!",xfile);

				else fclose(fpcheck);

				break;
		}
	}


	/*Checking SEISMOGRAM Output Particle velocities */
	/*-------------------------------------- */
	/* only MYID is performing the file checks, if to many PEs try to do it simultaneously,
	 * the code and/or FILESYSTEM and/or MPI implementation will cause segmentation faults
	 */
	if ((gv->SEISMO>0) && (MYID==0)) {

		fprintf(gv->FP," Check accessibility for seismogram files of each PE ... \n");
		fprintf(gv->FP," However, the list below refers only to PE0  ... \n");

		switch (gv->SEIS_FORMAT) {
			case 1:
				sprintf(file_ext,"su");
				break;

			case 2:
				sprintf(file_ext,"txt");
				break;

			case 3:
				sprintf(file_ext,"bin");
				break;
		}

		if (gv->SEIS_FORMAT==2) strcpy(xmod,"a");

		else strcpy(xmod,"w+b");

		/*MYID=0 is checking if all seismogram file written by other PEs can be written
		 * assuming that all PEs can address the files ystem equally
		 */
		for (myidcounter=0; myidcounter< (gv->NPROCX*gv->NPROCY); myidcounter++) {

			switch (gv->SEISMO) {
				case 1: /* particle velocities only */
					sprintf(xfile,"%s_vx.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					/*in case of number of PE's=500, there will be 500 messages, too many to be displayed! */
					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage);
					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					sprintf(xfile,"%s_vy.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					break;

				case 2 : /* pressure only */
					sprintf(xfile,"%s_p.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					break;

				case 4 : /* everything */
					sprintf(xfile,"%s_vx.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					/*in case of number of PE's=500, there will be 500 messages, too many to be displayed! */
					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage);
					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					sprintf(xfile,"%s_vy.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					//if (access(xfile,W_OK|X_OK)==-1) declare_error(errormessage);
					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					sprintf(xfile,"%s_p.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

				case 3 : /* curl and div only */
					sprintf(xfile,"%s_div.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					sprintf(xfile,"%s_curl.%s.%d",gv->SEIS_FILE,file_ext,myidcounter);

					if (MYID==myidcounter) fprintf(gv->FP,"    Check accessibility for seismogram file %s... \n",xfile);

					sprintf(errormessage,"PE %i cannot write seismogram file %s!",MYID,xfile);

					if ((fpcheck=fopen(xfile,xmod)) ==NULL) declare_error(errormessage);

					else fclose(fpcheck);

					remove(xfile);

					break;

			}
		}
	}

	/*Checking CHECKPOINT Output */
	/*-------------------------- */
	if (gv->CHECKPTREAD>0) {
		strcpy(xmod,"rb");
		sprintf(xfile,"%s.%d",gv->CHECKPTFILE,MYID);
		fprintf(gv->FP," Check readability for checkpoint files %s... \n",xfile);

		if (((fpcheck=fopen(xfile,xmod)) ==NULL) && (MYID==0)) declare_error(" PE 0 cannot read checkpoints!");

		else fclose(fpcheck);
	}

	if ((gv->CHECKPTWRITE>0)) {
		strcpy(xmod,"ab");
		sprintf(xfile,"%s.%d",gv->CHECKPTFILE,MYID);
		fprintf(gv->FP," Check writability for checkpoint files %s... \n",xfile);

		if (((fpcheck=fopen(xfile,xmod)) ==NULL) && (MYID==0)) declare_error(" PE 0 cannot write checkpoints!");

		else fclose(fpcheck);    /* Is there any reason to remove it? */
	}

	fprintf(gv->FP," Accessibility of output files from PE %d has been checked successfully.\n", MYID);

	fprintf(gv->FP,"\n\n --------- DETERMININATION OF MIN AND MAX VELOCITIES -----------\n");

	/* low Q frame not yet applied as a absorbing boundary */
	/* if (!FREE_SURF) ny1=1+FW;*/
	/*nfw=FW; check only outside the absorbing frame */
	nfw=0;
	cmax_s=0;
	cmin_s=10000;
	cmax_p=0;
	cmin_p=10000;


	/* find maximum model phase velocity of shear waves at infinite
	      frequency within the whole model */
	if (gv->L>0) {   /*viscoelastic*/
		for (i=1+nfw; i<= (nx-nfw); i++) {
			for (j=ny1; j<= (ny-nfw); j++) {
				c=sqrt(pu[j][i]* (1.0/prho[j][i]) * (1.0+gv->L*ptaus[j][i]));

				/*if c is close to zero (water, air), c will be ignored for finding
				cmax,cmin*/
				if ((cmax_s<c) && (c>cwater)) cmax_s=c;

				/* find minimum model phase velocity of shear waves at center
					       frequency of the source */
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

	/* find maximum model phase velocity of P-waves at infinite
		 frequency within the whole model */
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

	fprintf(gv->FP," Minimum and maximum P-wave and S-wave velocities within subvolumes: \n ");
	fprintf(gv->FP," MYID\t Vp_min(f=fc) \t Vp_max(f=inf) \t Vs_min(f=fc) \t Vs_max(f=inf) \n");
	fprintf(gv->FP," %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \n", MYID, cmin_p, cmax_p, cmin_s, cmax_s);

	fprintf(gv->FP," Note : if any P- or S-wave velocity is set below 1.0 m/s to simulate water or air,\n");
	fprintf(gv->FP," this minimum velocity will be ignored for determining stable DH and DT.\n\n");

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

	//if (gv->RSG) dtstab=gv->DH/cmax;

	if (MYID == 0) {

		fprintf(gv->FP," Global values for entire model: \n");
		fprintf(gv->FP," V_max= %5.2f m/s \t V_min= %5.2f m/s \n\n", cmax,cmin);
		fprintf(gv->FP,"\n\n ------------------ CHECK FOR GRID DISPERSION --------------------\n");
		fprintf(gv->FP," To satisfactorily limit grid dispersion the number of gridpoints \n");
		fprintf(gv->FP," per minimum wavelength (of S-waves) should be 6 (better more).\n");
		fprintf(gv->FP," Here the minimum wavelength is assumed to be minimum model phase velocity \n");
		fprintf(gv->FP," (of S-waves) at maximum frequency of the source\n");
		fprintf(gv->FP," devided by maximum frequency of the source.\n");
		fprintf(gv->FP," Maximum frequency of the source is approximately %6.3f Hz\n",2.0/gv->TS);
		fprintf(gv->FP," The minimum wavelength (P- or S-wave) in the following simulation will\n");
		fprintf(gv->FP," be %5.3f meter.\n", cmin/fmax);
		fprintf(gv->FP," Thus, the recommended value for DH is %5.3f meter.\n", dhstab);
		fprintf(gv->FP," You have specified DH= %5.3f meter.\n\n", gv->DH);

		if (gv->DH>dhstab)
			warning(" Grid dispersion will influence wave propagation, choose smaller grid spacing (DH).");


		fprintf(gv->FP," \n\n ----------------------- CHECK FOR STABILITY ---------------------\n");
		fprintf(gv->FP," The following simulation is stable provided that\n\n");

		//if (gv->RSG) fprintf(gv->FP," \t p=cmax*DT/DH < 1,\n\n");
		//else     fprintf(gv->FP," \t p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n\n");
        fprintf(gv->FP," \t p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n\n");

		fprintf(gv->FP," where cmax is the maximum phase velocity at infinite frequency\n");

		//if (gv->RSG==0) fprintf(gv->FP," and gamma = sum(|FD coeff.|)\n");
        fprintf(gv->FP," and gamma = sum(|FD coeff.|)\n");

		fprintf(gv->FP," In the current simulation cmax is %8.2f m/s .\n\n",cmax);

		fprintf(gv->FP," DT is the timestep and DH is the grid size.\n\n");
		fprintf(gv->FP," In this simulation the Courant-Friedrichs-Lewy number is %2.4f.\n",CFL);
		fprintf(gv->FP," In this simulation the stability limit for timestep DT is %e seconds .\n",dtstab);
		fprintf(gv->FP," You have specified DT= %e s.\n", gv->DT);

		if (gv->DT>dtstab)
			declare_error(" The simulation will get unstable, choose smaller DT. ");

		else fprintf(gv->FP," The simulation will be stable.\n");

		fprintf(gv->FP," \n\n --------------------- CHECK FOR INPUT ERRORS ---------------------\n");

		fprintf(gv->FP," Checking the time of wave propagation. \n");

		if (gv->SNAP) {
			fprintf(gv->FP," Checking the snapshot parameters. \n");

			if (gv->TSNAP2>gv->TIME) {
				sprintf(errormessage,"\nTSNAP2 = %e (last snapshot) > Time of wave propagation %e. TSNAP2 was changed to be equal to TIME.\n",gv->TSNAP2, gv->TIME);
				gv->TSNAP2=gv->TIME;

				if (MYID==0)
					warning(errormessage); /* if TSNAP2>simulation TIME, snapmerge will generate "additional" snapshots out of nowhere, thus, snapshot files size blow up */
			}

			snapoutx=gv->NX/ (float) gv->IDX;
			snapouty=gv->NY/ (float) gv->IDY;

			fprintf(gv->FP," Output of snapshot gridpoints per node (NX/NPROCX/IDX) %8.2f .\n", snapoutx);
			fprintf(gv->FP," Output of snapshot gridpoints per node (NY/NPROCY/IDY) %8.2f .\n", snapouty);

			if (snapoutx- (int) snapoutx>0)
				declare_error("\n\n Ratio NX-NPROCX-IDX must be whole-numbered \n\n");

			if (snapouty- (int) snapouty>0)
				declare_error("\n\n Ratio NY-NPROCY-IDY must be whole-numbered \n\n");

		}


		if ((gv->SEISMO) && (MYID==0)) {
			fprintf(gv->FP," Checking the number of seismogram samples. \n");
			fprintf(gv->FP,"    Number of timesteps %d.\n", gv->NT);
			fprintf(gv->FP,"    Seismogram sampling rate in timesteps %d.\n", gv->NDT);
			fprintf(gv->FP,"    Number of seismogram output samples %d.\n", gv->NT/gv->NDT);

			/* SU and SEG-Y allow 32767 samples, furthermore the exist programs allow for 65535 samples
			and pseudo-SEG-Y formats allow foralmost arbitrarily long traces.
			For binary and textual output the limit is arbitrary. USHRT_MAX is the limut of an unsigned short specified in limits.h */

			if ((gv->SEIS_FORMAT==1) && (gv->NT/gv->NDT) > (USHRT_MAX)) {
				fprintf(gv->FP," Maximum allowed number of samples per trace in SU format: %d \n", USHRT_MAX);
				declare_error(" Sorry. Too many samples per receiver! \n");
			}
		}

		if ((gv->SEISMO>0) && (MYID==0)) {
			srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
			srec_maxx=-1.0, srec_maxy=-1.0;
			fprintf(gv->FP," Checking for receiver position(s) specified in input file.\n");
			fprintf(gv->FP,"    Global grid size in m: %5.2f (x) : %5.2f (y) \n",gv->NX*gv->DH*gv->NPROCX,gv->NY*gv->DH*gv->NPROCY);

			if (gv->FREE_SURF==0) fprintf(gv->FP,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n", (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH, (float) gv->FW*gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);

			if (gv->FREE_SURF==1) fprintf(gv->FP,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n", (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH,gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);

			/* find maximum and minimum source positions coordinate ---- from input file*/
			/*for usability reasons, "z" - as commonly used - denotes the depth (vertical direction),
			      however, internally "y" is used for the vertical coordinate,
			      we simply switch the "y" and "z" coordinate as read in the input file,
			      therefore we determine the minimum/maximum position in y-direction by the ZREC1 variable and vice versa.
			this has to be considered for the receiver line coordinates specified in both the input file and separate source/receiver files*/

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

				fprintf(gv->FP,"    Number of receiver positions in input file : %i \n", ntr);
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

				fprintf(gv->FP,"    Number of receiver positions in receiver file %s : %i \n", gv->REC_FILE, ntr);
			}



			fprintf(gv->FP,"    Minimum receiver position coordinates : %5.2f (x) : %5.2f (y) \n",srec_minx,srec_miny);
			fprintf(gv->FP,"    Maximum receiver position coordinates : %5.2f (x) : %5.2f (y) \n",srec_maxx,srec_maxy);

			/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
			if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
				declare_error("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
			}

			if ((srec_maxx>gv->NX*gv->DH*gv->NPROCX) || (srec_maxy>gv->NY*gv->DH*gv->NPROCY)) {
				declare_error("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
			}

			/* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
			if ((srec_maxx< ((float) gv->FW*gv->DH)) || (srec_minx< ((float) gv->FW*gv->DH))) {
				/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
				warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 1, left boundary). \n\n");
			}

			if (srec_maxx> (gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 2, right boundary). \n\n");
			}

			if (srec_maxy> (gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 3, lower boundary). \n\n");
			}

			if ((srec_miny< ((float) gv->FW*gv->DH)) && !(gv->FREE_SURF)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 4, top boundary). \n\n");
			}

			fprintf(gv->FP," ... complete, receiver position specified in input file are located within the global grid.\n");

		}

		if ((gv->SRCREC==1) && (MYID==0)) {
			srec_minx=gv->DH*gv->NX*gv->NPROCX+1, srec_miny=gv->DH*gv->NY*gv->NPROCY+1;
			srec_maxx=-1.0, srec_maxy=-1.0;
			fprintf(gv->FP," Checking for source position(s) specified in source file. \n");
			fprintf(gv->FP,"    Global grid size in m: %5.2f (x) : %5.2f (y) \n",gv->NX*gv->DH*gv->NPROCX,gv->NY*gv->DH*gv->NPROCY);

			if (gv->FREE_SURF==0) fprintf(gv->FP,"    Global grid size in m (-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n", (float) gv->FW*gv->DH,gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH, (float) gv->FW*gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);

			if (gv->FREE_SURF==1) fprintf(gv->FP,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n", (float) gv->FW*gv->DH,gv->NX* (float) gv->DH*gv->NPROCX- (float) gv->FW*gv->DH,gv->DH,gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH);


			for (k=1; k<=nsrc; k++) {
				/* find maximum source positions coordinate*/
				if (srcpos[1][k]>srec_maxx) srec_maxx=srcpos[1][k];

				if (srcpos[2][k]>srec_maxy) srec_maxy=srcpos[2][k];

				/* find minimum source positions coordinate*/
				if (srcpos[1][k]<srec_minx) srec_minx=srcpos[1][k];

				if (srcpos[2][k]<srec_miny) srec_miny=srcpos[2][k];
			}

			fprintf(gv->FP,"    Number of source positions in source file %s. : %i.\n", gv->SOURCE_FILE, nsrc);
			fprintf(gv->FP,"    Minimum source position coordinates : %5.2f (x) : %5.2f (y) \n",srec_minx,srec_miny);
			fprintf(gv->FP,"    Maximum source position coordinates : %5.2f (x) : %5.2f (y) \n",srec_maxx,srec_maxy);

			/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
			if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
				declare_error("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
			}

			if ((srec_maxx>gv->NX*gv->DH*gv->NPROCX) || (srec_maxy>gv->NY*gv->DH*gv->NPROCY)) {
				declare_error("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
			}

			/* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
			if ((srec_maxx< ((float) gv->FW*gv->DH)) || (srec_minx< ((float) gv->FW*gv->DH))) {
				/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
				warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 1, left boundary). \n\n");
			}

			if (srec_maxx> (gv->NX*gv->DH*gv->NPROCX- (float) gv->FW*gv->DH)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 2, right boundary). \n\n");
			}

			if (srec_maxy> (gv->NY*gv->DH*gv->NPROCY- (float) gv->FW*gv->DH)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 3, lower boundary). \n\n");
			}

			if ((srec_miny< ((float) gv->FW*gv->DH)) && !(gv->FREE_SURF)) {
				/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
				warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 4, top boundary). \n\n");
			}

			fprintf(gv->FP," ...complete, all source position(s) specified in source file are located within the global grid.\n");


		}

		fprintf(gv->FP,"\n\n ----------------------- ABSORBING BOUNDARY ------------------------\n");



		if (gv->ABS_TYPE==1) {
			fprintf(gv->FP," You have specified a CPML boundary (ABS_TYPE=1) with a width of %i gridpoints (%5.2f m).\n",gv->FW, (float) gv->FW*gv->DH);

			if (gv->FW<10) {
				fprintf(gv->FP," Width (FW) of absorbing frame should be at least 10 gridpoints.\n\n");
				warning("  Be aware of artificial reflections from grid boundaries ! \n");

			}
		}

		if (gv->ABS_TYPE==2) {
			fprintf(gv->FP," You have specified an absorbing boundary (ABS_TYPE=2) width of %i gridpoints (%5.2f m).\n",gv->FW, (float) gv->FW*gv->DH);

			if (gv->FW<30) {
				fprintf(gv->FP," Width (FW) of absorbing frame should be at least 30 gridpoints.\n\n");
				warning(" Be aware of artificial reflections from grid boundaries ! \n");
			}
		}

		if (((gv->NX) < gv->FW) || ((gv->NY) <gv->FW))	{
			fprintf(gv->FP," \n Width of absorbing boundary (FW = %i gridpoints) is larger than at least one subdomain dimension: \n",gv->FW);
			fprintf(gv->FP," \t NX/NPROCX = %i, NY/NPROCY = %i gridpoints.\n",gv->NX,gv->NY);
			declare_error(" Choose smaller width of absorbing frame (FW) or increase subdomain dimensions");
		}

		if (gv->BOUNDARY) {
			if (gv->ABS_TYPE==1 || gv->ABS_TYPE==2) {
				warning(" You have activated a periodic boundary and set an absorbing boundary at the same time! \n");
			}
        }
	}
}
