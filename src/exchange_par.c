
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
 * Exchange FD-Parameters between PEs
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "enums.h"

void exchange_par(GlobVar *gv, GlobVarInv *vinv)
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

        fdum[42] = vinv->EPSILON;
        fdum[43] = vinv->EPSILON_u;
        fdum[44] = vinv->EPSILON_rho;
        fdum[45] = vinv->SRTRADIUS;
        fdum[46] = vinv->VPUPPERLIM;
        fdum[47] = vinv->VPLOWERLIM;
        fdum[48] = vinv->VSUPPERLIM;
        fdum[49] = vinv->VSLOWERLIM;
        fdum[50] = vinv->RHOUPPERLIM;
        fdum[51] = vinv->RHOLOWERLIM;
        fdum[52] = vinv->F_LOW_PASS_START;
        fdum[53] = vinv->F_LOW_PASS_END;
        fdum[54] = vinv->F_LOW_PASS_INCR;
        fdum[55] = vinv->EPS_SCALE;
        fdum[56] = vinv->SCALEFAC;
        fdum[57] = vinv->PRO;
        fdum[58] = vinv->TWLENGTH_PLUS;
        fdum[59] = vinv->TWLENGTH_MINUS;
        fdum[60] = vinv->GAMMA;
        //
        //
        fdum[63] = vinv->WATERLEVEL_LNORM8;
        fdum[64] = vinv->VP_VS_RATIO;
        fdum[65] = vinv->S_VS;
        fdum[66] = vinv->S_VP;
        fdum[67] = vinv->S_RHO;
        fdum[68] = vinv->A;
        fdum[69] = vinv->F_HIGH_PASS;
        //fdum[70] = vinv->JOINT_INVERSION_PSV_SH_ALPHA_VS;
        //fdum[71] = vinv->JOINT_INVERSION_PSV_SH_ALPHA_RHO;
        fdum[72] = vinv->EPSILON_WE;
        // fdum[73] = vinv->EPSILON_WE_SH;
        fdum[74] = vinv->WOLFE_C1_SL;
        fdum[75] = vinv->WOLFE_C2_SL;
        fdum[76] = vinv->TRKILL_STF_OFFSET_LOWER;
        fdum[77] = vinv->TRKILL_STF_OFFSET_UPPER;
        fdum[78] = vinv->TRKILL_OFFSET_LOWER;
        fdum[79] = vinv->TRKILL_OFFSET_UPPER;
        fdum[80] = vinv->LBFGS_SCALE_GRADIENTS;

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
        idum[10] = gv->SOURCE_TOPO;
        idum[11] = gv->READMOD;
        idum[12] = gv->L;
        idum[13] = gv->FREE_SURF;
        idum[14] = gv->SNAP;
        idum[15] = gv->SNAPSHOT_START;
        idum[16] = gv->SNAPSHOT_END;
        idum[17] = gv->SNAPSHOT_INCR;
        idum[18] = gv->DRX;
        idum[19] = gv->BOUNDARY;
        idum[20] = gv->REC_ARRAY;
        idum[21] = gv->SRCREC;
        idum[22] = gv->IDX;
        idum[23] = gv->IDY;
        idum[24] = gv->MODE;
        idum[25] = gv->WEQ;
        idum[26] = gv->SNAP_FORMAT;
        idum[27] = gv->SEISMO;
        idum[28] = gv->READREC;
        idum[29] = gv->REC_TOPO;
        idum[30] = gv->NDT;
        idum[31] = gv->SEIS_FORMAT;
        idum[32] = gv->NT;
        idum[33] = gv->NS;
        idum[34] = gv->FDORDER;
        idum[35] = gv->MAXRELERROR;
        idum[36] = gv->RUN_MULTIPLE_SHOTS;
        idum[37] = gv->WRITE_MODELFILES;
        idum[38] = gv->ABS_TYPE;
        idum[39] = gv->FDORDER_TIME;
        idum[40] = gv->SIGOUT;
        idum[41] = gv->SIGOUT_FORMAT;
        idum[42] = gv->SOURCE_SHAPE_OLD;

/* ****** FWI ****** */
        idum[45] = vinv->TRKILL_STF;
        idum[46] = vinv->GRADT1;
        idum[47] = vinv->GRADT2;
        idum[48] = vinv->GRADT3;
        idum[49] = vinv->GRADT4;
        idum[50] = vinv->ITERMAX;
        //
        //
        idum[54] = vinv->ADJOINT_TYPE;
        idum[55] = vinv->TESTSHOT_START;
        idum[56] = vinv->TESTSHOT_END;
        idum[57] = vinv->TESTSHOT_INCR;
        idum[58] = vinv->SWS_TAPER_GRAD_VERT;
        idum[59] = vinv->SWS_TAPER_GRAD_HOR;
        idum[60] = vinv->SWS_TAPER_GRAD_SOURCES;
        idum[61] = vinv->SWS_TAPER_CIRCULAR_PER_SHOT;
        idum[62] = vinv->SRTSHAPE;
        idum[63] = vinv->FILTSIZE;
        idum[64] = vinv->SPATFILTER;
        idum[65] = vinv->SPAT_FILT_SIZE;
        idum[66] = vinv->SPAT_FILT_1;
        idum[67] = vinv->SPAT_FILT_ITER;
        idum[68] = vinv->INV_RHO_ITER;
        idum[69] = vinv->NFSTART;
        idum[70] = vinv->NF;
        idum[71] = vinv->NFSTART_JAC;
        idum[72] = vinv->NF_JAC;
        idum[73] = vinv->SWS_TAPER_FILE;
        idum[75] = vinv->GRAD_METHOD;
        idum[76] = vinv->MODEL_FILTER;
        idum[77] = vinv->FILT_SIZE;
        //
        idum[79] = vinv->INV_STF;
        idum[80] = vinv->N_STF;
        idum[81] = vinv->N_STF_START;
        idum[82] = vinv->TIME_FILT;
        idum[83] = vinv->ORDER;
        idum[84] = vinv->LNORM;
        idum[85] = vinv->DTINV;
        idum[86] = vinv->NTDTINV;
        idum[87] = vinv->STEPMAX;
        idum[88] = vinv->TRKILL;
        idum[89] = vinv->TIMEWIN;
        idum[90] = vinv->NORMALIZE;
        idum[91] = vinv->INV_VP_ITER;
        idum[92] = vinv->INV_VS_ITER;
        idum[93] = vinv->MIN_ITER;
        idum[94] = vinv->GRAD_FILTER;
        idum[95] = vinv->FILT_SIZE_GRAD;
        idum[96] = vinv->NO_OF_TESTSHOTS;
        idum[97] = vinv->VELOCITY;
        idum[98] = vinv->SWS_TAPER_FILE_PER_SHOT;
        idum[99] = vinv->S;
        idum[100] = vinv->GRAD_FILT_WAVELENGTH;
        idum[101] = vinv->WRITE_DIFF;
        //
        //  idum[103] = vinv->WAVETYPE;
        //  idum[104] = vinv->SOURCE_SHAPE_SH;
        //  idum[105] = vinv->JOINT_INVERSION_PSV_SH_TYPE;
        //
        //
        //
        //
        idum[110] = vinv->USE_WORKFLOW;
        idum[111] = vinv->EPRECOND;
        idum[112] = vinv->EPRECOND_ITER;
        //
        idum[114] = vinv->LBFGS_STEP_LENGTH;
        idum[115] = vinv->EPRECOND_PER_SHOT;
        //  idum[116] = vinv->EPRECOND_PER_SHOT_SH;
        idum[117] = vinv->N_LBFGS;
        idum[118] = vinv->WOLFE_CONDITION;
        idum[119] = vinv->WOLFE_NUM_TEST;
        idum[120] = vinv->WOLFE_TRY_OLD_STEPLENGTH;
        idum[121] = vinv->WRITE_FILTERED_DATA;
        idum[122] = vinv->TAPER_STF;
        idum[123] = vinv->TW_IND;
        idum[124] = vinv->TRKILL_OFFSET;
        idum[125] = vinv->TRKILL_STF_OFFSET;
        idum[126] = vinv->TRKILL_STF_OFFSET_INVERT;
        //  idum[127] = vinv->JOINT_EQUAL_WEIGHTING;
        idum[128] = vinv->STF_FULL;
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

    MPI_Bcast(&(vinv->DATA_DIR), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->MISFIT_LOG_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->INV_MODELFILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->JACOBIAN), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->FILE_WORKFLOW), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->WORKFLOW_HEADER), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->PARA), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->TRKILL_FILE_STF), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->TAPER_FILE_NAME), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->FREQ_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(vinv->TRKILL_FILE), STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

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

/* ****** FWI ****** */

    vinv->EPSILON = fdum[42];
    vinv->EPSILON_u = fdum[43];
    vinv->EPSILON_rho = fdum[44];
    vinv->SRTRADIUS = fdum[45];
    vinv->VPUPPERLIM = fdum[46];
    vinv->VPLOWERLIM = fdum[47];
    vinv->VSUPPERLIM = fdum[48];
    vinv->VSLOWERLIM = fdum[49];
    vinv->RHOUPPERLIM = fdum[50];
    vinv->RHOLOWERLIM = fdum[51];
    vinv->F_LOW_PASS_START = fdum[52];
    vinv->F_LOW_PASS_END = fdum[53];
    vinv->F_LOW_PASS_INCR = fdum[54];
    vinv->EPS_SCALE = fdum[55];
    vinv->SCALEFAC = fdum[56];
    vinv->PRO = fdum[57];
    vinv->TWLENGTH_PLUS = fdum[58];
    vinv->TWLENGTH_MINUS = fdum[59];
    vinv->GAMMA = fdum[60];
    //
    //
    vinv->WATERLEVEL_LNORM8 = fdum[63];
    vinv->VP_VS_RATIO = fdum[64];
    vinv->S_VS = fdum[65];
    vinv->S_VP = fdum[66];
    vinv->S_RHO = fdum[67];
    vinv->A = fdum[68];
    vinv->F_HIGH_PASS = fdum[69];
    //vinv->JOINT_INVERSION_PSV_SH_ALPHA_VS = fdum[70];
    //vinv->JOINT_INVERSION_PSV_SH_ALPHA_RHO = fdum[71];
    vinv->EPSILON_WE = fdum[72];
    //vinv->EPSILON_WE_SH = fdum[73];
    vinv->WOLFE_C1_SL = fdum[74];
    vinv->WOLFE_C2_SL = fdum[75];
    vinv->TRKILL_STF_OFFSET_LOWER = fdum[76];
    vinv->TRKILL_STF_OFFSET_UPPER = fdum[77];
    vinv->TRKILL_OFFSET_LOWER = fdum[78];
    vinv->TRKILL_OFFSET_UPPER = fdum[79];
    vinv->LBFGS_SCALE_GRADIENTS = fdum[80];

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
    gv->SOURCE_TOPO = idum[10];
    gv->READMOD = idum[11];
    gv->L = idum[12];
    gv->FREE_SURF = idum[13];
    gv->SNAP = idum[14];
    gv->SNAPSHOT_START = idum[15];
    gv->SNAPSHOT_END = idum[16];
    gv->SNAPSHOT_INCR = idum[17];
    gv->DRX = idum[18];
    gv->BOUNDARY = idum[19];
    gv->REC_ARRAY = idum[20];
    gv->SRCREC = idum[21];
    gv->IDX = idum[22];
    gv->IDY = idum[23];
    gv->MODE = idum[24];
    gv->WEQ = (WEQTYPE)idum[25];
    gv->SNAP_FORMAT = idum[26];
    gv->SEISMO = idum[27];
    gv->READREC = idum[28];
    gv->REC_TOPO = idum[29];
    gv->NDT = idum[30];
    gv->SEIS_FORMAT = idum[31];
    gv->NT = idum[32];
    gv->NS = idum[33];
    gv->FDORDER = idum[34];
    gv->MAXRELERROR = idum[35];
    gv->RUN_MULTIPLE_SHOTS = idum[36];
    gv->WRITE_MODELFILES = idum[37];
    gv->ABS_TYPE = idum[38];
    gv->FDORDER_TIME = idum[39];
    gv->SIGOUT = idum[40];
    gv->SIGOUT_FORMAT = idum[41];
    gv->SOURCE_SHAPE_OLD = idum[42];

/* ****** FWI ****** */

    vinv->TRKILL_STF = idum[45];
    vinv->GRADT1 = idum[46];
    vinv->GRADT2 = idum[47];
    vinv->GRADT3 = idum[48];
    vinv->GRADT4 = idum[49];
    vinv->ITERMAX = idum[50];
    //
    //  
    vinv->ADJOINT_TYPE = idum[54];
    vinv->TESTSHOT_START = idum[55];
    vinv->TESTSHOT_END = idum[56];
    vinv->TESTSHOT_INCR = idum[57];
    vinv->SWS_TAPER_GRAD_VERT = idum[58];
    vinv->SWS_TAPER_GRAD_HOR = idum[59];
    vinv->SWS_TAPER_GRAD_SOURCES = idum[60];
    vinv->SWS_TAPER_CIRCULAR_PER_SHOT = idum[61];
    vinv->SRTSHAPE = idum[62];
    vinv->FILTSIZE = idum[63];
    vinv->SPATFILTER = idum[64];
    vinv->SPAT_FILT_SIZE = idum[65];
    vinv->SPAT_FILT_1 = idum[66];
    vinv->SPAT_FILT_ITER = idum[67];
    vinv->INV_RHO_ITER = idum[68];
    vinv->NFSTART = idum[69];
    vinv->NF = idum[70];
    vinv->NFSTART_JAC = idum[71];
    vinv->NF_JAC = idum[72];
    vinv->SWS_TAPER_FILE = idum[73];
    vinv->GRAD_METHOD = idum[75];
    vinv->MODEL_FILTER = idum[76];
    vinv->FILT_SIZE = idum[77];
    //  
    vinv->INV_STF = idum[79];
    vinv->N_STF = idum[80];
    vinv->N_STF_START = idum[81];
    vinv->TIME_FILT = idum[82];
    vinv->ORDER = idum[83];
    vinv->LNORM = idum[84];
    vinv->DTINV = idum[85];
    vinv->NTDTINV = idum[86];
    vinv->STEPMAX = idum[87];
    vinv->TRKILL = idum[88];
    vinv->TIMEWIN = idum[89];
    vinv->NORMALIZE = idum[90];
    vinv->INV_VP_ITER = idum[91];
    vinv->INV_VS_ITER = idum[92];
    vinv->MIN_ITER = idum[93];
    vinv->GRAD_FILTER = idum[94];
    vinv->FILT_SIZE_GRAD = idum[95];
    vinv->NO_OF_TESTSHOTS = idum[96];
    vinv->VELOCITY = idum[97];
    vinv->SWS_TAPER_FILE_PER_SHOT = idum[98];
    vinv->S = idum[99];
    vinv->GRAD_FILT_WAVELENGTH = idum[100];
    vinv->WRITE_DIFF = idum[101];
    //
    //    vinv->WAVETYPE = idum[103];
    //    vinv->SOURCE_SHAPE_SH = idum[104];
    //    vinv->JOINT_INVERSION_PSV_SH_TYPE = idum[105];
    //
    //
    //
    //  
    vinv->USE_WORKFLOW = idum[110];
    vinv->EPRECOND = idum[111];
    vinv->EPRECOND_ITER = idum[112];
    //  
    vinv->LBFGS_STEP_LENGTH = idum[114];
    vinv->EPRECOND_PER_SHOT = idum[115];
    //  vinv->EPRECOND_PER_SHOT_SH = idum[116];
    vinv->N_LBFGS = idum[117];
    vinv->WOLFE_CONDITION = idum[118];
    vinv->WOLFE_NUM_TEST = idum[119];
    vinv->WOLFE_TRY_OLD_STEPLENGTH = idum[120];
    vinv->WRITE_FILTERED_DATA = idum[121];
    vinv->TAPER_STF = idum[122];
    vinv->TW_IND = idum[123];
    vinv->TRKILL_OFFSET = idum[124];
    vinv->TRKILL_STF_OFFSET = idum[125];
    vinv->TRKILL_STF_OFFSET_INVERT = idum[126];
    //  vinv->JOINT_EQUAL_WEIGHTING = idum[127];
    vinv->STF_FULL = idum[128];

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
