/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------
 * Merge and output of the gradient
 *-------------------------------------------------------------------------------*/

#include "fd.h"

void output_gradient(int ishot, int nshots, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    char modfile[2*STRING_SIZE];
    
    /* write gradient file - vp*/
    sprintf(modfile,"%s_gradient.vp.shot%i",vinv->GRADIENT,ishot);
    writemod(modfile, minv->gradVp_shot, 3, gv);
    MPI_Barrier(MPI_COMM_WORLD);
    if (gv->MPID==0) mergemod(modfile, 3, gv);

    /* write gradient file - vs*/
    sprintf(modfile,"%s_gradient.vs.shot%i",vinv->GRADIENT,ishot);
    writemod(modfile, minv->gradVs_shot, 3, gv);
	MPI_Barrier(MPI_COMM_WORLD);
	if (gv->MPID==0) mergemod(modfile, 3, gv);

    /* write gradient file - rho'*/
	sprintf(modfile,"%s_gradient.rho.shot%i",vinv->GRADIENT,ishot);
    writemod(modfile, minv->gradRho_shot, 3, gv);
    MPI_Barrier(MPI_COMM_WORLD);
    if (gv->MPID==0) mergemod(modfile, 3, gv);

    /* output of gradient summed up for all shots */
    if (ishot == nshots){

        /* write gradient file - vp*/
        sprintf(modfile,"%s_gradient.vp",vinv->GRADIENT);
        writemod(modfile, minv->gradVp, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID==0) mergemod(modfile, 3, gv);

        /* write gradient file - vs*/
        sprintf(modfile,"%s_gradient.vs",vinv->GRADIENT);
        writemod(modfile, minv->gradVs, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID==0) mergemod(modfile, 3, gv);

        /* write gradient file - rho'*/
        sprintf(modfile,"%s_gradient.rho",vinv->GRADIENT);
        writemod(modfile, minv->gradRho, 3, gv);
        MPI_Barrier(MPI_COMM_WORLD);
        if (gv->MPID==0) mergemod(modfile, 3, gv);
    }
}
