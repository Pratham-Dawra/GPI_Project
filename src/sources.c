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
/* ----------------------------------------------------------------------
 * Reading (distributed) source positions, timeshift, centre frequency
 * and amplitude from SOURCE_FILE.
 *
 * ---------------------------------------------------------------------- */

#include "fd.h"

float **sources(int *nsrc, GlobVar *gv){

	/* declaration of extern variables */
	//extern float PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE, TS, DH;
	//extern char SOURCE_FILE[STRING_SIZE];
	//extern int MYID, NXG, NYG, SRCREC, RUN_MULTIPLE_SHOTS, SOURCE_TYPE, SOURCE_SHAPE;
	//extern FILE *FP;

	float **srcpos=NULL;
	int   i, l, isrc=0, current_source=0,nvarin=0;
	float xsrc, ysrc, tshift, tan_phi, dz, x, fc=0.0;
	char buffer[STRING_SIZE], bufferstring[10], cline[256];
	FILE *fpsrc;

    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

	if (MYID==0){
		fprintf(gv->FP," Message from function sources (written by PE %d):\n",MYID);
		if (gv->SRCREC==1){ /* read source positions from file */
			*nsrc=0;
			fprintf(gv->FP,"\n ------------------ READING SOURCE PARAMETERS ------------------- \n");
			fprintf(gv->FP,"\n Reading source parameters from file: %s (SOFI2D source format)\n",gv->SOURCE_FILE);
			if ((fpsrc=fopen(gv->SOURCE_FILE,"r"))==NULL) declare_error(" Source file could not be opened !");
			while(fgets(buffer, STRING_SIZE, fpsrc))
			{
				sscanf(buffer,"%s",bufferstring);
				/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
				/*if (sscanf(buffer,"%s",bufferstring)==1) printf("string %s \n",bufferstring);*/
				if ((strchr(buffer,'#')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(*nsrc);
			}

			rewind(fpsrc);

			if ((nsrc)==0) fprintf(gv->FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %i.\n",(*nsrc=0));
			else fprintf(gv->FP," Number of source positions specified in %s : %d \n",gv->SOURCE_FILE,*nsrc);

			srcpos=matrix(1,NSPAR,1,*nsrc);

			/* memory for source position definition (Ricker, Fuchs-Mueller, sin**3 & from_File) */
			if (gv->SOURCE_SHAPE <= 4) {
				/* srcpos[1][l] = x position
				   srcpos[2][l] = depth position
				   srcpos[3][l] = horizontal position (always zero in 2D)
				   srcpos[4][l] = time delay (source start time)
				   srcpos[5][l] = centre frequency
				   srcpos[6][l] = amplitude
				   srcpos[7][l] = azimuth [°] (optional)
				   srcpos[8][l] = SOURCE_TYPE (optional)
				*/

				for (l=1;l<=*nsrc;l++){
					fgets(cline,255,fpsrc);
					nvarin=sscanf(cline,"%f%f%f%f%f%f%f",&xsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l], &srcpos[8][l]);
					switch(nvarin){
					case 0: xsrc=0.0;
					case 1: ysrc=0.0;
					case 2: if (MYID==0) fprintf(gv->FP," No time shift defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 3: if (MYID==0) fprintf(gv->FP," No frequency defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 4: if (MYID==0) fprintf(gv->FP," No amplitude defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 5: srcpos[7][l]=0.0;
					case 6: srcpos[8][l]=gv->SOURCE_TYPE;
					}
					if ((srcpos[8][l]!=4) && (nvarin>5)) {
						current_source=(int)srcpos[8][l];
						if (MYID==0) fprintf(gv->FP," SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n", l, current_source);
					}
					/* fscanf(fpsrc,"%f%f%f%f%f",&xsrc, &ysrc, &tshift, &fc, &amp); */

					srcpos[1][l]=xsrc;
					srcpos[2][l]=ysrc;
					srcpos[3][l]=0.0;
					srcpos[4][l]=tshift;
					fc=srcpos[5][l];
				}
			}
			/* memory for source position definition (Berlage) */
			else if (5 == gv->SOURCE_SHAPE) {
				/* srcpos[1][l] = x position
				   srcpos[2][l] = depth position
				   srcpos[3][l] = horizontal position (always zero in 2D)
				   srcpos[4][l] = time delay (source start time)
				   srcpos[5][l] = centre frequency
				   srcpos[6][l] = amplitude
				   srcpos[7][l] = azimuth [°] (optional)
				   srcpos[8][l] = SOURCE_TYPE (optional)
				   srcpos[9][l] = time exponent (Berlage only)
				   srcpos[10][l] = exponential decay factor (Berlage only)
				   srcpos[11][l] = initial phase angle [°] (Berlage only)
				*/

				for (l=1;l<=*nsrc;l++){
					fgets(cline,255,fpsrc);
					nvarin=sscanf(cline,"%f%f%f%f%f%f%f%f%f%f",&xsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[9][l], &srcpos[10][l], &srcpos[11][l], &srcpos[7][l], &srcpos[8][l]);
					switch(nvarin){
					case 0: xsrc=0.0;
					case 1: ysrc=0.0;
					case 2: if (MYID==0) fprintf(gv->FP," No time shift defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 3: if (MYID==0) fprintf(gv->FP," No frequency defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 4: if (MYID==0) fprintf(gv->FP," No amplitude defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 5: if (MYID==0) fprintf(gv->FP," No time exponent defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 6: if (MYID==0) fprintf(gv->FP," No exponential decay factor defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 7: if (MYID==0) fprintf(gv->FP," No initial phase angle [°] defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 8: srcpos[7][l]=0.0;
					case 9: srcpos[8][l]=gv->SOURCE_TYPE;
					}
				if ((srcpos[8][l]!=4) && (nvarin>8)) {
					current_source=(int)srcpos[8][l];
					if (MYID==0) fprintf(gv->FP," SOURCE_TYPE of source #%i is specified as %i,    SOURCE_AZIMUTH is ignored.\n", l, current_source);
				}

				srcpos[1][l]=xsrc;
				srcpos[2][l]=ysrc;
				srcpos[3][l]=0.0;
				srcpos[4][l]=tshift;
				fc=srcpos[5][l];
				}
			}
			/* memory for source position definition (Klauder) */
			else if (6 == gv->SOURCE_SHAPE) {
				/* srcpos[1][l] = x position
				   srcpos[2][l] = depth position
				   srcpos[3][l] = horizontal position (always zero in 2D)
				   srcpos[4][l] = time delay (source start time)
				   srcpos[5][l] = centre frequency
				   srcpos[6][l] = amplitude
				   srcpos[7][l] = azimuth [°] (optional)
				   srcpos[8][l] = SOURCE_TYPE (optional)
				   srcpos[9][l] = minimum frequency (Klauder only)
				   srcpos[10][l] = maximum frequency (Klauder only)
				   srcpos[11][l] = sweep length [s] (Klauder only)
				   srcpos[12][l] = width of the Klauder wavelet in number of centre periods (Klauder only)
				*/

				for (l=1;l<=*nsrc;l++){
					fgets(cline,255,fpsrc);
					nvarin=sscanf(cline,"%f%f%f%f%f%f%f%f%f%f",&xsrc, &ysrc, &tshift, &srcpos[9][l], &srcpos[10][l], &srcpos[6][l], &srcpos[11][l], &srcpos[12][l], &srcpos[7][l], &srcpos[8][l]);
					switch(nvarin){
					case 0: xsrc=0.0;
					case 1: ysrc=0.0;
					case 2: if (MYID==0) fprintf(gv->FP," No time shift defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 3: if (MYID==0) fprintf(gv->FP," No minimum frequency defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 4: if (MYID==0) fprintf(gv->FP," No maximum frequency defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 5: if (MYID==0) fprintf(gv->FP," No amplitude defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 6: if (MYID==0) fprintf(gv->FP," No sweep length [s] defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 7: if (MYID==0) fprintf(gv->FP," No width of the wavelet (in number of centre periods) defined for source %i in %s!\n",l, gv->SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 8: srcpos[7][l]=0.0;
					case 9: srcpos[8][l]=gv->SOURCE_TYPE;
					}
					if ((srcpos[8][l]!=4) && (nvarin>8)) {
						current_source=(int)srcpos[8][l];
						if (MYID==0) fprintf(gv->FP," SOURCE_TYPE of source #%i is specified as %i,    SOURCE_AZIMUTH is ignored.\n", l, current_source);
					}

					srcpos[1][l]=xsrc;
					srcpos[2][l]=ysrc;
					srcpos[3][l]=0.0;
					srcpos[4][l]=tshift;
					fc=(srcpos[9][l]+srcpos[10][l])/2;
					srcpos[5][l] = fc;
				}
			}

			fclose(fpsrc);

			/* Compute maximum frequency */
			for (l=1;l<=*nsrc;l++)
				if (srcpos[5][l]>fc) fc=srcpos[5][l];
			fprintf(gv->FP," Maximum frequency defined in %s: %6.2f Hz\n",gv->SOURCE_FILE,fc);
			gv->TS=1.0/fc;

			/* outputs all sources per each subdomain / node*/

			if (MYID==0){
				/*fprintf(gv->FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");

				for (l=1;l<=*nsrc;l++)
					fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n\n",
							l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);*/
				if (gv->RUN_MULTIPLE_SHOTS) fprintf(gv->FP," All sources will be modelled individually because of RUN_MULTIPLE_SHOTS=1!\n\n");
				else fprintf(gv->FP," All sources will be modelled simultaneously because of RUN_MULTIPLE_SHOTS=0!\n\n");

			}

		} 
		else if (gv->SRCREC==2) {
			if (gv->PLANE_WAVE_DEPTH > 0) {  /* plane wave excitation */
				if (gv->SOURCE_SHAPE > 4) {
					declare_error("Plane wave is only implemented for Ricker, Fuchs-Mueller, sin**3 or an external wavelet! Change parameter SOURCE_SHAPE!");
				}

				fprintf(gv->FP," Computing source nodes for plane wave excitation.\n");
				fprintf(gv->FP," depth= %5.2f meter, incidence angle= %5.2f degrees.\n",gv->PLANE_WAVE_DEPTH, gv->PLANE_WAVE_ANGLE);


				tan_phi=tan(gv->PLANE_WAVE_ANGLE*PI/180.0);

				dz=(float)gv->NXG*gv->DH*tan_phi;
				fprintf(gv->FP," Message from function sources (written by PE %d):\n",MYID);
				fprintf(gv->FP," Maximum depth of plane wave: %5.2f meter \n",gv->PLANE_WAVE_DEPTH+dz);
				if ((gv->PLANE_WAVE_DEPTH+dz)<=gv->NYG*gv->DH){
					*nsrc=gv->NXG;
					srcpos=matrix(1,8,1,*nsrc);
					isrc=0;
					for (i=1;i<=gv->NXG;i++){
						isrc++;
						x=(float)i*gv->DH;
						srcpos[1][isrc]=x;
						srcpos[2][isrc]=gv->PLANE_WAVE_DEPTH+(tan_phi*x);
						srcpos[3][isrc]=0.0;
						srcpos[4][isrc]=0.0;
						srcpos[5][isrc]=1.0/gv->TS;
						srcpos[6][isrc]=1.0;
						srcpos[7][isrc]=0.0;
						srcpos[8][isrc]=gv->SOURCE_TYPE;
					}
				}
				else declare_error(" Maximum depth of plane wave exceeds model depth. ");
			}
			else declare_error("SRCREC parameter specifies PLANE_WAVE excitation, but PLANE_WAVE_DEPTH<=0!");
		}
		else declare_error("SRCREC parameter is invalid (SRCREC!=1 or SRCREC!=2)! No source parameters specified!");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&(gv->TS),1,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (MYID!=0) srcpos=matrix(1,NSPAR,1,*nsrc);
	MPI_Bcast(&srcpos[1][1],(*nsrc)*12,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (MYID==0){
		if (*nsrc>50) fprintf(gv->FP," The following table is quite large (%i lines) and will, thus, be truncated to the first 50 entries! \n\n",*nsrc);
		if (4 <= gv->SOURCE_SHAPE) {
			fprintf(gv->FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");
		}
		else if (5 == gv->SOURCE_SHAPE) {
			fprintf(gv->FP," number\t  x\t\t  y\t\t  tshift\t  fc\t\t  amp\t  n\t  alpha\t  phi0\t  source_azimuth\t  source_type\n");
		}
		else if (6 == gv->SOURCE_SHAPE) {
			fprintf(gv->FP," number\t  x\t\t  y\t\t  tshift\t  fmin\t  fmax\t  fc\t\t  amp\t  dt_sweep\t  dt_wavelet\t  source_azimuth\t  source_type\n");
		}

		if (*nsrc>50) { for (l=1;l<=50;l++)
			if (4 <= gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);
			}
			else if (5 == gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[9][l],srcpos[10][l],srcpos[11][l],srcpos[7][l],srcpos[8][l]);
			}
			else if (6 == gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[9][l],srcpos[10][l],srcpos[5][l],srcpos[6][l],srcpos[11][l],srcpos[12][l],srcpos[7][l],srcpos[8][l]);
			}
		}
		else for (l=1;l<=*nsrc;l++) {
			if (4 <= gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);
			}
			else if (5 == gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[9][l],srcpos[10][l],srcpos[11][l],srcpos[7][l],srcpos[8][l]);
			}
			else if (6 == gv->SOURCE_SHAPE) {
				fprintf(gv->FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %1.0f\n", l,srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[9][l],srcpos[10][l],srcpos[5][l],srcpos[6][l],srcpos[11][l],srcpos[12][l],srcpos[7][l],srcpos[8][l]);
			}
		}
	}

	return srcpos;
}
