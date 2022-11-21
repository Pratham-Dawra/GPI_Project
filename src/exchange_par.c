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
 *   Exchange FD-Parameters between PEs                         
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(GlobVar *gv) {

  /* definition of local variables */
  int idum[NPAR];
  float fdum[NPAR];
    
  int MYID;
  MPI_Comm_rank(MPI_COMM_WORLD, &MYID);
  
  if (MYID == 0) {
    fdum[1]  = gv->DH;
    fdum[2]  = gv->TIME;
    fdum[3]  = gv->DT;
    fdum[4]  = gv->TS;
    fdum[5]  = 0.0;
    fdum[6]  = 0.0;

    fdum[8]  = gv->TAU;

    fdum[10]  = gv->TSNAP1;
    fdum[11]  = gv->TSNAP2;
    fdum[12]  = gv->TSNAPINC;
    fdum[13]  = gv->REFREC[1];
    fdum[14]  = gv->REFREC[2];
    fdum[15]  = gv->PLANE_WAVE_ANGLE;

    fdum[16]  = gv->XREC1;
    fdum[17]  = gv->YREC1;

    fdum[19]  = gv->XREC2;
    fdum[20]  = gv->YREC2;

    fdum[22]  = gv->DAMPING;
    fdum[23]  = gv->REC_ARRAY_DEPTH;
    fdum[24]  = gv->REC_ARRAY_DIST;
    fdum[25]  = gv->PLANE_WAVE_DEPTH;

    fdum[26]  = gv->NGEOPH;

    fdum[27]  = 0.0; //gv->SRCPOSXYZ[0];
    fdum[28]  = 0.0; //gv->SRCPOSXYZ[1];
    fdum[29]  = 0.0; //gv->SRCPOSXYZ[2];

    fdum[30]  = gv->FPML;
    fdum[31]  = gv->VPPML;
    fdum[32]  = gv->NPOWER;
    fdum[33]  = gv->K_MAX_CPML;

    /*************************************/

    idum[1]  = gv->NPROCX;
    idum[2]  = gv->NPROCY;
    idum[3]  = gv->LOG;

    idum[4]  = gv->NPROC;
    idum[5]  = gv->NX;
    idum[6]  = gv->NY;
    idum[7]  = gv->FW;
    idum[8]  = gv->SOURCE_SHAPE;
    idum[9]  = gv->SOURCE_TYPE;
    idum[10]  = gv->READMOD;
    idum[11]  = gv->L;
    idum[12]  = gv->FREE_SURF;
    idum[13]  = gv->SNAP;
    idum[14]  = gv->DRX;

    idum[16]  = gv->BOUNDARY;
    idum[17]  = gv->REC_ARRAY;
    idum[18]  = gv->SRCREC;
    idum[19]  = gv->IDX;
    idum[20]  = gv->IDY;
    idum[21]  = 0;
    idum[22]  = gv->WEQ;
    idum[23]  = gv->SNAP_FORMAT;
    idum[24]  = gv->SEISMO;
    idum[25]  = gv->READREC;
    idum[26]  = 0; //gv->RSG;
    idum[27]  = gv->NDT;
    idum[28]  = gv->SEIS_FORMAT;
    idum[29]  = gv->CHECKPTREAD;
    idum[30]  = gv->CHECKPTWRITE;

    idum[31]  = gv->FDORDER;
    idum[32]  = gv->MAXRELERROR;
    idum[33]  = gv->RUN_MULTIPLE_SHOTS;
    idum[34]  = gv->WRITE_MODELFILES;

    idum[35] = gv->ABS_TYPE;

    idum[36] = gv->FDORDER_TIME;

  } /** if (MYID == 0) **/

  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&(gv->SOURCE_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->MFILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->SNAP_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->REC_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->SEIS_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->LOG_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->SIGNAL_FILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&(gv->CHECKPTFILE),STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  gv->DH = fdum[1];
  gv->TIME = fdum[2];
  gv->DT = fdum[3];
  gv->TS = fdum[4];

  gv->TAU = fdum[8];

  gv->TSNAP1 = fdum[10];
  gv->TSNAP2 = fdum[11];
  gv->TSNAPINC = fdum[12];
  gv->REFREC[1] = fdum[13];
  gv->REFREC[2] = fdum[14];
  gv->PLANE_WAVE_ANGLE = fdum[15];
  gv->XREC1 = fdum[16];
  gv->YREC1 = fdum[17];

  gv->XREC2 = fdum[19];
  gv->YREC2 = fdum[20];

  gv->DAMPING = fdum[22];
  gv->REC_ARRAY_DEPTH = fdum[23];
  gv->REC_ARRAY_DIST = fdum[24];
  gv->PLANE_WAVE_DEPTH = fdum[25];

  gv->NGEOPH = fdum[26];

  //SRCPOSXYZ[0] = fdum[27];
  //SRCPOSXYZ[1] = fdum[28];
  //SRCPOSXYZ[2] = fdum[29];

  gv->FPML = fdum[30];
  gv->VPPML = fdum[31];
  gv->NPOWER = fdum[32];
  gv->K_MAX_CPML = fdum[33];

  /********************************************/

  gv->NPROCX = idum[1];
  gv->NPROCY = idum[2];
  gv->LOG=idum[3];
  gv->NPROC  = idum[4];
  gv->NX = idum[5];
  gv->NY = idum[6];
  gv->FW = idum [7];
  gv->SOURCE_SHAPE = idum[8];
  gv->SOURCE_TYPE = idum[9];
  gv->READMOD = idum[10];
  gv->L = idum[11];
  gv->FREE_SURF = idum[12];
  gv->SNAP = idum[13];
  gv->DRX = idum[14];

  gv->BOUNDARY = idum[16];
  gv->REC_ARRAY = idum[17];
  gv->SRCREC = idum[18];
  gv->IDX = idum[19];
  gv->IDY = idum[20];
  gv->WEQ= idum[22];

  gv->SNAP_FORMAT = idum[23];
  gv->SEISMO = idum[24];
  gv->READREC = idum[25];
  //RSG = idum[26];
  gv->NDT = idum[27];
  gv->SEIS_FORMAT = idum[28];
  gv->CHECKPTREAD = idum[29];
  gv->CHECKPTWRITE = idum[30];

  gv->FDORDER = idum[31];
  gv->MAXRELERROR = idum[32];
  gv->RUN_MULTIPLE_SHOTS = idum[33];
  gv->WRITE_MODELFILES = idum[34];

  gv->ABS_TYPE = idum[35];

  gv->FDORDER_TIME = idum[36];
  
  if (gv->L>0) {
    if (MYID != 0) {
      gv->FL = vector(1,gv->L);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&gv->FL[1],gv->L,MPI_FLOAT,0,MPI_COMM_WORLD);
  } else {
    gv->FL = NULL;
  }
}
