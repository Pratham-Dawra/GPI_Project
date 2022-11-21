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
 *   Write FD-Parameters to stdout                           
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp, GlobVar *gv) {

	/* definition of local variables */
	int l;
	char file_ext[8];
    int MYID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
    
	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ********* PARAMETERS AS SPECIFIED IN INPUT FILE **********");
	fprintf(fp,"\n **********************************************************\n\n");

	fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
	fprintf(fp,"\n");
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",gv->NPROCX);
	fprintf(fp," Number of PEs in vertical direction (NPROCY): %d\n",gv->NPROCY);
	fprintf(fp," Total number of PEs in use: %d\n",gv->NP);
	fprintf(fp,"\n");

	fprintf(fp,"------------------------- FD Algorithm ------------------------\n");
	/*if (gv->RSG) {
		fprintf(fp," Rotated Staggered Grid (RSG) is used. \n");
		fprintf(fp," Order of spatial FD operators (FDORDER) is set to 2.");

	} else 	{
		fprintf(fp," Standard Staggered Grid (SSG) (Virieux-grid) is used.\n");
		fprintf(fp," Order of spatial FD operators (FDORDER) is %d\n", gv->FDORDER);
        fprintf(fp," Order of temporal FD operator (FDORDER_TIME) is %d\n", gv->FDORDER_TIME);
	}*/
    fprintf(fp," Standard Staggered Grid (SSG) (Virieux-grid) is used.\n");
    fprintf(fp," Order of spatial FD operators (FDORDER) is %d\n", gv->FDORDER);
    fprintf(fp," Order of temporal FD operator (FDORDER_TIME) is %d\n", gv->FDORDER_TIME);
	fprintf(fp,"\n");

    fprintf(fp,"----------------------- Wave Equation ------------------------\n");

    switch (gv->WEQ){
        case 1 :
            fprintf(fp," Acoustic wave equation\n");
            break;
        case 2 :
            fprintf(fp," Viscoacoustic wave equation\n");
            break;
        case 3 :
            fprintf(fp," Elastic wave equation\n");
            break;
        case 4 :
            fprintf(fp," Viscoelastic wave equation\n");
            break;
        case 5 :
            fprintf(fp," Elastic VTI wave equation\n");
            break;
        case 6 :
            fprintf(fp," Viscoelastic VTI wave equation\n");
            break;
        case 7 :
            fprintf(fp," Elastic TTI wave equation\n");
            break;
        case 8 :
            fprintf(fp," Viscoelastic TTI wave equation\n");
            break;
            
        default :
            declare_error(" Sorry, incorrect specification of type of wave equation (WEQ) ! ");
            break;
    }

    fprintf(fp,"\n");
    

	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", gv->NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", gv->NY);
	fprintf(fp," Grid-spacing (DH): %e meter\n", gv->DH);
	fprintf(fp," Model size in x-direction: %.5g m\n", gv->NX*gv->DH);
	fprintf(fp," Model size in y-direction: %.5g m\n", gv->NY*gv->DH);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",gv->TIME);
	fprintf(fp," Timestep (DT): %e seconds\n", gv->DT);
	fprintf(fp," Number of timesteps: %i \n",gv->NT);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");

	if (gv->SRCREC) {
		fprintf(fp," reading source positions, time delay, centre frequency \n");
		fprintf(fp," and initial amplitude from ASCII-file: %s\n\n",gv->SOURCE_FILE);
	} else {
		fprintf(fp," plane wave excitation: depth= %5.2f meter \n",gv->PLANE_WAVE_DEPTH);
		fprintf(fp," incidence angle of plane P-wave (from vertical) PLANE_WAVE_ANGLE= %5.2f degrees \n",gv->PLANE_WAVE_ANGLE);
 		fprintf(fp," duration of source signal: %e seconds\n",gv->TS);
 		fprintf(fp," (centre frequency is approximately %e Hz)\n",1.0/gv->TS);
	}

	fprintf(fp," wavelet of source:");

	switch (gv->SOURCE_SHAPE) {
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," reading from \n\t %s\n",gv->SIGNAL_FILE);
		break;
	case 4 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		break;
	case 5 :
		fprintf(fp," Berlage wavelet \n");
		break;
	case 6 :
		fprintf(fp," Klauder wavelet \n");
		break;
	default :
		declare_error(" Sorry, incorrect specification of source wavelet ! ");
		break;
	}

	fprintf(fp," Type of source:");
	switch (gv->SOURCE_TYPE) {
	case 1 :
		fprintf(fp," explosive source \n");
		break;
	case 2 :
		fprintf(fp," point source with directive force in x-direction\n");
		break;
	case 3 :
		fprintf(fp," point source with directive force in (vertical) y-direction\n");
		break;
	case 4 :
		fprintf(fp," point source with directive force in  z-direction\n");
		break;
	default :
		declare_error(" Sorry, wrong source type specification ! ");
		break;
	}

	if (1==gv->SIGOUT) {
        switch (gv->SIGOUT_FORMAT){
            case 1: sprintf(file_ext,"su");  break;
            case 2: sprintf(file_ext,"txt"); break;
            case 3: sprintf(file_ext,"bin"); break;
		}
        fprintf(fp," Source signal will be written to file: %s.%s\n\n",gv->SIGOUT_FILE,file_ext);
    }

	fprintf(fp,"\n");

	if (gv->SEISMO){
		fprintf(fp," ------------------------- RECEIVER  ------- -------------------\n");
		if (gv->READREC){
			fprintf(fp," reading receiver positions from file \n");
			fprintf(fp,"\t%s\n\n",gv->REC_FILE);
			fprintf(fp," reference point for receiver coordinate system:\n");
			fprintf(fp," x=%f \ty=%f\t z=%f\n",gv->REFREC[1], gv->REFREC[2], gv->REFREC[3]);
		} else if (gv->REC_ARRAY>0){
				fprintf(fp," horizontal lines of receivers.\n");
				fprintf(fp," number of lines: %d \n",gv->REC_ARRAY);
				fprintf(fp," depth of upper line: %e m \n",gv->REC_ARRAY_DEPTH);
				fprintf(fp," vertical increment between lines: %e m \n",gv->REC_ARRAY_DIST);
				fprintf(fp," distance between receivers in x-direction within line: %i \n", gv->DRX);		
		}else{

			fprintf(fp," first receiver position (XREC1,YREC1) = (%e, %e) m\n",
			    gv->XREC1,gv->YREC1);
			fprintf(fp," last receiver position (XREC2,YREC2) = (%e, %e) m\n",
			    gv->XREC2,gv->YREC2);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
	if (gv->FREE_SURF) fprintf(fp," free surface at the top of the model ! \n");
	else fprintf(fp," no free surface at the top of the model ! \n");
	fprintf(fp,"\n");

	fprintf(fp," ------------------------- ABSORBING FRAME ---------------------\n");
	if (gv->FW>0){
		fprintf(fp," width of absorbing frame is %i grid points (%5.2f m).\n",gv->FW,(float)gv->FW*gv->DH);
		if (gv->ABS_TYPE==1) {
                	fprintf(fp," CPML damping applied. \n");
                	fprintf(fp," Damping velocity in the PML frame in m/s: %f .\n",gv->VPPML);
                	fprintf(fp," Frequency within the PML frame in Hz: %f \n",gv->FPML);
                	fprintf(fp," NPOWER: %f \n",gv->NPOWER);
                	fprintf(fp," K_MAX: %f \n",gv->K_MAX_CPML);
        	}

		if (gv->ABS_TYPE==2) {
			fprintf(fp," Exponential damping applied. \n");
			fprintf(fp," Percentage of amplitude decay: %f .\n",gv->DAMPING);
		}
	}
	else fprintf(fp," absorbing frame not installed ! \n");

	switch (gv->BOUNDARY){
	case 0 :
		fprintf(fp," No periodic boundary condition.\n");
		break;
	case 1 :
		fprintf(fp," Periodic boundary condition at left and right edges.\n");
		break;
	default :
		warning(" Wrong integer value for BOUNDARY specified ! ");
		warning(" No periodic boundary condition will be applied ");
		gv->BOUNDARY=0;
		break;
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
	fprintf(fp," Number of relaxation mechanisms (L): %i\n",gv->L);
	fprintf(fp," The L relaxation frequencies are at:  \n");
	for (l=1;l<=gv->L;l++) fprintf(fp,"\t%f",gv->FL[l]);
	fprintf(fp," Hz\n");
	fprintf(fp," Value for tau is : %f\n",gv->TAU);


	if (gv->SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(gv->SNAP){
		case 1:
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			break;
		case 2:
			fprintf(fp," pressure field.\n");
			break;
		case 3:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			break;
		case 4:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			fprintf(fp," x- and y-component of particle velocity.\n");
			break;
		default:
			declare_error(" sorry, incorrect value for SNAP ! \n");
			break;
		}

		fprintf(fp," \t first (TSNAP1)= %8.5f s\n", gv->TSNAP1);
		fprintf(fp," \t last (TSNAP2)=%8.5f s\n",gv->TSNAP2);
		fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",gv->TSNAPINC);
		fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,gv->NX);
		fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,gv->NY);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",gv->SNAP_FILE);
		switch (gv->SNAP_FORMAT){
		case 1 :
			declare_error(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			declare_error(" Don't know the format for the Snapshot-data ! \n");
			break;
		}

		fprintf(fp,"\n\n");
	}
	if (gv->SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		switch (gv->SEIS_FORMAT){
			case 1: sprintf(file_ext,"su");  break;
			case 2: sprintf(file_ext,"txt"); break;
			case 3: sprintf(file_ext,"bin"); break;
		}
		if ((gv->SEISMO==1) || (gv->SEISMO==4)){
			fprintf(fp," Seismograms of ");
			fprintf(fp," x-, y-, and z-component");
			fprintf(fp," of particle velocity.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_vx.%s\n\t%s_vy.%s\n",gv->SEIS_FILE,file_ext,gv->SEIS_FILE,file_ext);
		}
		if ((gv->SEISMO==2) || (gv->SEISMO==4)){
			fprintf(fp," Seismograms of pressure field (hydrophones).\n");
			fprintf(fp," output-file: \n ");
			fprintf(fp,"\t%s_p.%s\n",gv->SEIS_FILE,file_ext);
		}
		if ((gv->SEISMO==3) || (gv->SEISMO==4)){
			fprintf(fp," Seismograms of curl (S-wave component) and div (P-wave component of wavefield).\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_rot.%s \n\t%s_div.%s\n",gv->SEIS_FILE,file_ext,gv->SEIS_FILE,file_ext);
			
		}		

		switch (gv->SEIS_FORMAT){
		case 1 :
			fprintf(fp," The data is written in IEEE SU-format . \n");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary IEEE (4 byte per float)");
			break;
		default:
			declare_error(" Sorry. I don't know the format for the seismic data ! \n");
			break;
		}
		fprintf(fp," samplingrate of seismic data: %f s\n",gv->NDT*gv->DT);
		if (!gv->READREC) fprintf(fp," Trace-spacing: %5.2f m\n", gv->NGEOPH*gv->DH);
		fprintf(fp," Number of samples per trace: %i \n", iround(gv->NT/gv->NDT));
		fprintf(fp," ----------------------------------------------------------\n");
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ******* PARAMETERS READ or PROCESSED within SOFI2D ********");
	fprintf(fp,"\n **********************************************************\n\n");
}
