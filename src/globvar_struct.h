
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
 * globvar_struct.h - global variables of viscoelastic 2D FD program
 * generally, for the names of the global variables uppercase letters are used
 *
 * ----------------------------------------------------------------------*/

#ifndef GLOBVAR_STRUCT_H_INCLUDED
#define GLOBVAR_STRUCT_H_INCLUDED

#include "macros.h"

typedef void (*FDop_s_fct)(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy);
typedef void (*FDop_v_fct)(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy, float **sxy);

typedef struct {
//struct ModelVar {
    // Models
    int READMOD;                // switch to read model parameters from MFILE
    char MFILE[STRING_SIZE];    // model file name
    int NX;                     // number of grid points in x-direction
    int NY;                     // number of grid points in y-direction (depth)
    /*Attenuation */
    int L;                      // number of relaxation parameters
    float *FL;                  // frequency of each relaxation parameters [Hz]
    float TAU;                  // ratio of retardation and relaxation time
    /*Source */
    int SOURCE_TYPE;            // type of source
    int SOURCE_SHAPE;           // shape of source-signal
    char SIGNAL_FILE[STRING_SIZE];  // name of external signal file
    int SRCREC;                 // switch to read source parameters from external source file
    char SOURCE_FILE[STRING_SIZE];  // name of source parameter file
    int RUN_MULTIPLE_SHOTS;     // multiple shots modeled simultaneously (0) or individually; added for multiple shots
    float PLANE_WAVE_DEPTH;     // depth of plane wave excitation [meter]
    float PLANE_WAVE_ANGLE;     // dip of plane wave from vertical [degree]
    float TS;                   // duration of source signal [in second]
    float XS;                   // Source location ??
    float YS;                   // Source location ??
    /*Receiver */
    int READREC;                // switch to read receiver positions from file
    char REC_FILE[STRING_SIZE]; // name of external receiver file
    float REFREC[4];            // reference point for receiver coordinate system
    float XREC1;                // x-position of first receiver [m]
    float XREC2;                // x-position of last receiver [m]
    float YREC1;                // y-position of first receiver [m]
    float YREC2;                // y-position of last receiver [m]
    float NGEOPH;               // distance between two adjacent receivers [gridpoints]; in auto mode NGEOPH will be 
    // calculated from model dimensions, type integer is incorrect
    int REC_ARRAY;              // number of receivers in 1D receiver array
    float REC_ARRAY_DEPTH;      // depth of first plane [m] 
    float REC_ARRAY_DIST;       // increment between receiver planes [m]
    int DRX;                    // increment between receivers in each plane [gridpoints]
//}

//struct IOVar {
    //IO
    /*Seismograms */
    int SEISMO;                 // switch to output components of seismograms
    int NDT;                    // sampling rate of seismograms [timesteps DT]
    int SEIS_FORMAT;            // data output format for seismograms
    char SEIS_FILE[STRING_SIZE];    // name of output file of seismograms
    /*Snapshots */
    int SNAP;                   // switch to output of snapshots
    int SNAP_FORMAT;            // data output format for snapshots
    char SNAP_FILE[STRING_SIZE];    // name of output file of snapshots
    float TSNAP1;               // first snapshot [s]
    float TSNAP2;               // last snapshot [s]
    float TSNAPINC;             // increment between snapshots [s]
    int IDX;                    // increment in x-direction [gridpoints]
    int IDY;                    // increment in y-direction [gridpoints]
    /*Others */
    int WRITE_MODELFILES;       // switch to output model files
    int SIGOUT;                 // switch to output source wavelet
    int SIGOUT_FORMAT;          // data output format for source wavelet
    char SIGOUT_FILE[STRING_SIZE];  // name of output file of source wavelet
    int LOG;                    // switch to output logging information
    char LOG_VERBOSITY[STRING_SIZE];    // log output level (verbosity)       
    char LOG_FILE[STRING_SIZE]; // name of output file of logging information
    int OUTNTIMESTEPINFO;       // every OUTNTIMESTEPINFO th timestep, information on the time step will be given to screen/file
    int CHECKPTREAD;            // switch to read wavefield from checkpoint file
    int CHECKPTWRITE;           // switch to save wavefield to checkpoint file
    char CHECKPTFILE[STRING_SIZE];  // name of checkpoint file
//}

//struct FDParams{
    // FD Params
    int WEQ;                    // wave equation
    float DH;                   // spacial increment [m]
    int FDORDER;                // spatial FD order
    float TIME;                 // time (of modelling) [s]
    float DT;                   // time increment (of modelling) [s]
    int FDORDER_TIME;           // temporal FD order
    int NT;                     //// number of timesteps (=iround(TIME/DT))
    FDop_s_fct FDOP_S;          // function pointer for FD operator
    FDop_v_fct FDOP_V;          // function pointer for FD operator
    /* MPI-variables */
    int MAXRELERROR;            // switch of maximum relative group velocity error
    int NPROCX;                 // number of processors in x-direction
    int NPROCY;                 // number of processors in y-direction
    int NPROC;                  //// number of processors (=NPROCX*NPROCY)
    int NP;                     //// number of processors from mpirun command
    int MPID;                   //// ID of processor
    int INDEX[5];               //// ID of neighboring processes %% Why 5??? - should be just 4
    int POS[3];                 //// processor location in the 3D logical processor array (MPID%NPROCX; MPID/NPROCX) %% Why 3?? - should be just 2
    int IENDX;                  //// size of domain in x-direction (=NX/NPROCX)
    int IENDY;                  //// size of domain in y-direction (=NY/NPROCY)
    int NXG;                    //// number of grid points in x-direction (global)
    int NYG;                    //// number of grid points in y-direction (global)
    const int TAG1;             //// ??? place holder ???
    const int TAG2;             //// ??? place holder ???
    const int TAG3;             //// ??? place holder ???
    const int TAG4;             //// ??? place holder ???
    const int TAG5;             //// ??? place holder ???
    const int TAG6;             //// ??? place holder ???
    // Boundary
    int FREE_SURF;              // switch to apply free surface at the top of the model
    int BOUNDARY;               // switch to apply periodic boundary condition at edges
    int ABS_TYPE;               // type of the absorbing boundary
    int FW;                     // width of absorbing frame [gridpoints]
    float DAMPING;              // attenuation at the edges of the grid [%]
    /* PML-Parameters */
    float NPOWER;               // exponent for calculation of damping profile
    float K_MAX_CPML;           // 
    float FPML;                 // dominant signal frequency (usually FC) [Hz]
    float VPPML;                // attenuation velocity within the PML boundary [m/s]
} GlobVar;

/* Parameters not used - please double check! */
//    int   RSG;                          //// ??? NOT USED ??? - rotated staggered grid; spatial adaptive Code variables
//    float SRCPOSXYZ[3]={0.0, 0.0, 0.0}; //// !!! NOT USED !!!
//    int   NSRC;                         //// !!! NOT USED !!!
//    int   NPSP;                         //// !!! NOT USED !!!
//    int   DC;                           //// !!! NOT USED !!!
//    int   check_id;                     //// !!! NOT USED !!!
//    int   cfgt_id;                      //// !!! NOT USED !!!
//    int   cfgt;                         //// !!! NOT USED !!!
//    int   jumpid;                       //// !!! NOT USED !!!
//    float DH1;                          //// !!! NOT USED !!!

#endif
