
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

/* ------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "enums.h"

void exchange_par(GlobVar * gv)
{
    int idum[NPAR];
    float fdum[NPAR];

    if (gv->MPID == 0) {
        fdum[1] = gv->DH;
        fdum[2] = gv->TIME;
        fdum[3] = gv->DT;
        fdum[4] = gv->TS;
	//
	//
	//
        fdum[8] = gv->TAU;
        fdum[9] = gv->F_REF;
        fdum[10] = gv->TSNAP1;
        fdum[11] = gv->TSNAP2;
        fdum[12] = gv->TSNAPINC;
        fdum[13] = gv->REFREC[1];
        fdum[14] = gv->REFREC[2];
        fdum[15] = gv->PLANE_WAVE_ANGLE;
        fdum[16] = gv->XREC1;
        fdum[17] = gv->YREC1;
	//
        fdum[19] = gv->XREC2;
        fdum[20] = gv->YREC2;
	//
        fdum[22] = gv->DAMPING;
        fdum[23] = gv->REC_ARRAY_DEPTH;
        fdum[24] = gv->REC_ARRAY_DIST;
        fdum[25] = gv->PLANE_WAVE_DEPTH;
        fdum[26] = gv->NGEOPH;
	//
	//
	//
        fdum[30] = gv->FPML;
        fdum[31] = gv->VPPML;
        fdum[32] = gv->NPOWER;
        fdum[33] = gv->K_MAX_CPML;

/* ****** FWI ****** */       

        // fdum[26]  = MUN;
        fdum[42]  = EPSILON;
        fdum[43]  = EPSILON_u;
        fdum[44]  = EPSILON_rho;
        
        fdum[45]  = SRTRADIUS;
        
        fdum[46]  = VPUPPERLIM;
        fdum[47]  = VPLOWERLIM;
        fdum[48]  = VSUPPERLIM;
        fdum[49]  = VSLOWERLIM;
        fdum[50]  = RHOUPPERLIM;
        fdum[51]  = RHOLOWERLIM;

        fdum[52]  = F_LOW_PASS_START;
        fdum[53]  = F_LOW_PASS_END;
        fdum[54]  = F_LOW_PASS_INCR;
        
        fdum[55]  = EPS_SCALE;
        fdum[56]  = SCALEFAC;
        fdum[57]  = PRO;
        
        fdum[58]  = TWLENGTH_PLUS;
        fdum[59]  = TWLENGTH_MINUS;
        fdum[60]  = GAMMA;
        
        fdum[62]  = TSHIFT_back;
        
        fdum[63]  = WATERLEVEL_LNORM8;
        
        fdum[64]  = VP_VS_RATIO;
        fdum[65]  = S_VS;
        fdum[66]  = S_VP;
        fdum[67]  = S_RHO;
        fdum[68]  = A;
        ´
        fdum[69]  = F_HIGH_PASS;
        
        fdum[70] = JOINT_INVERSION_PSV_SH_ALPHA_VS;
        fdum[71] = JOINT_INVERSION_PSV_SH_ALPHA_RHO;
        
        fdum[72]=EPSILON_WE;
        fdum[73]=EPSILON_WE_SH;
        
        fdum[74]=WOLFE_C1_SL;
        fdum[75]=WOLFE_C2_SL;
        
        fdum[76]=TRKILL_STF_OFFSET_LOWER;
        fdum[77]=TRKILL_STF_OFFSET_UPPER;
        fdum[78]=TRKILL_OFFSET_LOWER;
        fdum[79]=TRKILL_OFFSET_UPPER;
        
        fdum[80]=LBFGS_SCALE_GRADIENTS;

	/*************************************/

        idum[1] = gv->NPROCX;
        idum[2] = gv->NPROCY;
        idum[3] = gv->LOG;
        idum[4] = gv->NPROC;
        idum[5] = gv->NXG;
        idum[6] = gv->NYG;
        idum[7] = gv->FW;
        idum[8] = gv->SOURCE_SHAPE;
        idum[9] = gv->SOURCE_TYPE;
        idum[10] = gv->READMOD;
        idum[11] = gv->L;
        idum[12] = gv->FREE_SURF;
        idum[13] = gv->SNAP;
        idum[14] = gv->DRX;
	//
        idum[16] = gv->BOUNDARY;
        idum[17] = gv->REC_ARRAY;
        idum[18] = gv->SRCREC;
        idum[19] = gv->IDX;
        idum[20] = gv->IDY;
	//
        idum[22] = gv->WEQ;
        idum[23] = gv->SNAP_FORMAT;
        idum[24] = gv->SEISMO;
        idum[25] = gv->READREC;
	//
        idum[27] = gv->NDT;
        idum[28] = gv->SEIS_FORMAT;
        idum[29] = gv->NT;
        idum[30] = gv->NS;
        idum[31] = gv->FDORDER;
        idum[32] = gv->MAXRELERROR;
        idum[33] = gv->RUN_MULTIPLE_SHOTS;
        idum[34] = gv->WRITE_MODELFILES;
        idum[35] = gv->ABS_TYPE;
        idum[36] = gv->FDORDER_TIME;
        idum[37] = gv->SIGOUT;
        idum[38] = gv->SIGOUT_FORMAT;
	//

/* ****** FWI ****** */
        idum[21]  = TRKILL_STF;
        idum[34]  = TAPERLENGTH;
        idum[35]  = INVTYPE;
        idum[36]  = GRADT1;
        idum[37]  = GRADT2;
        idum[38]  = GRADT3;
        idum[39]  = GRADT4;
        idum[40]  = ITERMAX;
        idum[41]  = PARAMETERIZATION;
        idum[42]  = FW;
        idum[43]  = FORWARD_ONLY;
        idum[44]  = ADJOINT_TYPE;
        
        idum[45]  = TESTSHOT_START;
        idum[46]  = TESTSHOT_END;
        idum[47]  = TESTSHOT_INCR;
        
        idum[48]  = SWS_TAPER_GRAD_VERT;
        idum[49]  = SWS_TAPER_GRAD_HOR;
        idum[50]  = SWS_TAPER_GRAD_SOURCES;
        idum[51]  = SWS_TAPER_CIRCULAR_PER_SHOT;
        idum[52]  = SRTSHAPE;
        idum[53]  = FILTSIZE;
        
        idum[54]  = SPATFILTER;
        idum[55]  = SPAT_FILT_SIZE;
        idum[56]  = SPAT_FILT_1;
        idum[57]  = SPAT_FILT_ITER;
        
        idum[58]  = INV_RHO_ITER;
        idum[59]  = nfstart;
        idum[60]  = nf;
        
        idum[61]  = nfstart_jac;
        idum[62]  = nf_jac;
        idum[63]  = SWS_TAPER_FILE;
        idum[65]  = GRAD_METHOD;
        
        idum[66]  = MODEL_FILTER;
        idum[67]  = FILT_SIZE;
        
        idum[69]  = INV_STF;
        idum[70]  = N_STF;
        idum[71]  = N_STF_START;
        
        idum[72]  = TIME_FILT;
        idum[73]  = ORDER;
        
        idum[74]  = LNORM;
        idum[75]  = DTINV;
        
        idum[76]  = STEPMAX;
        
        idum[77]  = TRKILL;
        
        idum[78]  = TIMEWIN;
        
        idum[79]  = NORMALIZE;
        
        idum[80]  = INV_VP_ITER;
        idum[81]  = INV_VS_ITER;
        
        idum[82]  = MIN_ITER;
        
        idum[83]  = GRAD_FILTER;
        idum[84]  = FILT_SIZE_GRAD;
        
        idum[85]  = NO_OF_TESTSHOTS;
        
        // idum[86]  = EMPTY;
        
        idum[87]  = VELOCITY;
        
        idum[88]  = SWS_TAPER_FILE_PER_SHOT;
        
        idum[89]  = S;
        
        idum[90]  = GRAD_FILT_WAVELENGTH;
        
        idum[91]  = ACOUSTIC;
        
        
        idum[92]  = VERBOSE;
        
        idum[93]  = WAVETYPE;
        
        idum[94]  = SOURCE_SHAPE_SH;
        
        idum[95] = JOINT_INVERSION_PSV_SH_TYPE;
        
        idum[96]  = SNAPSHOT_START;
        idum[97]  = SNAPSHOT_END;
        idum[98]  = SNAPSHOT_INCR;
        
        idum[100]=USE_WORKFLOW;
        
        idum[101]=EPRECOND;
        idum[102]=EPRECOND_ITER;
        
        idum[104]=LBFGS_STEP_LENGTH;
        
        idum[105]=EPRECOND_PER_SHOT;
        idum[106]=EPRECOND_PER_SHOT_SH;
        
        idum[107]=N_LBFGS;

        idum[108]=WOLFE_CONDITION;
        idum[109]=WOLFE_NUM_TEST;
        idum[110]=WOLFE_TRY_OLD_STEPLENGTH;
        
        idum[111]=WRITE_FILTERED_DATA;
        
        idum[112]=TAPER_STF;
        idum[113]=TW_IND;
        
        idum[114]=TRKILL_OFFSET;
        idum[115]=TRKILL_STF_OFFSET;
        idum[116]=TRKILL_STF_OFFSET_INVERT;
        
        idum[117]=JOINT_EQUAL_WEIGHTING;
        idum[118]=STF_FULL;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&idum, NPAR, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fdum, NPAR, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&(gv->SOURCE_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->MFILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->SNAP_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->REC_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->SEIS_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->LOG_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->SIGNAL_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->SIGOUT_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->LOG_VERBOSITY), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    gv->DH = fdum[1];
    gv->TIME = fdum[2];
    gv->DT = fdum[3];
    gv->TS = fdum[4];
    //
    //
    //
    gv->TAU = fdum[8];
    gv->F_REF = fdum[9];
    gv->TSNAP1 = fdum[10];
    gv->TSNAP2 = fdum[11];
    gv->TSNAPINC = fdum[12];
    gv->REFREC[1] = fdum[13];
    gv->REFREC[2] = fdum[14];
    gv->PLANE_WAVE_ANGLE = fdum[15];
    gv->XREC1 = fdum[16];
    gv->YREC1 = fdum[17];
    // 
    gv->XREC2 = fdum[19];
    gv->YREC2 = fdum[20];
    //
    gv->DAMPING = fdum[22];
    gv->REC_ARRAY_DEPTH = fdum[23];
    gv->REC_ARRAY_DIST = fdum[24];
    gv->PLANE_WAVE_DEPTH = fdum[25];
    gv->NGEOPH = fdum[26];
    //
    //
    //
    gv->FPML = fdum[30];
    gv->VPPML = fdum[31];
    gv->NPOWER = fdum[32];
    gv->K_MAX_CPML = fdum[33];

    /********************************************/

    gv->NPROCX = idum[1];
    gv->NPROCY = idum[2];
    gv->LOG = idum[3];
    gv->NPROC = idum[4];
    gv->NXG = idum[5];
    gv->NYG = idum[6];
    gv->FW = idum[7];
    gv->SOURCE_SHAPE = idum[8];
    gv->SOURCE_TYPE = idum[9];
    gv->READMOD = idum[10];
    gv->L = idum[11];
    gv->FREE_SURF = idum[12];
    gv->SNAP = idum[13];
    gv->DRX = idum[14];
    //
    gv->BOUNDARY = idum[16];
    gv->REC_ARRAY = idum[17];
    gv->SRCREC = idum[18];
    gv->IDX = idum[19];
    gv->IDY = idum[20];
    //
    gv->WEQ = (WEQTYPE)idum[22];
    gv->SNAP_FORMAT = idum[23];
    gv->SEISMO = idum[24];
    gv->READREC = idum[25];
    //
    gv->NDT = idum[27];
    gv->SEIS_FORMAT = idum[28];
    gv->NT = idum[29];
    gv->NS = idum[30];
    gv->FDORDER = idum[31];
    gv->MAXRELERROR = idum[32];
    gv->RUN_MULTIPLE_SHOTS = idum[33];
    gv->WRITE_MODELFILES = idum[34];
    gv->ABS_TYPE = idum[35];
    gv->FDORDER_TIME = idum[36];
    gv->SIGOUT = idum[37];
    gv->SIGOUT_FORMAT = idum[38];

    if (gv->L > 0) {
        if (gv->MPID != 0) {
            gv->FL = vector(1, gv->L);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&gv->FL[1], gv->L, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else {
        gv->FL = NULL;
    }
}
