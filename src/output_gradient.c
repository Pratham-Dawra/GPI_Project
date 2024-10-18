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

    int i, j;
    char modfile[2*STRING_SIZE];

    FILE *FP6 = NULL;


    /* prepare gradient output: merge gradient files */
    /* merge gradient file - vp */
    sprintf(modfile,"%s_gradient_vp_shot.%i.%i.%i",vinv->GRADIENT,ishot,gv->POS[1],gv->POS[2]);
	FP6=fopen(modfile,"wb");

	for (i=1;i<=gv->NX;i=i+gv->IDX){
		for (j=1;j<=gv->NY;j=j+gv->IDY){
			fwrite(&minv->gradVp_shot[j][i],sizeof(float),1,FP6);
		}
	}

	fclose(FP6);

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(modfile,"%s_gradient_vp_shot.%i",vinv->GRADIENT,ishot);
	if (gv->MPID==0) mergemod(modfile, 3, gv);

    writemod(modfile, minv->gradVp_shot, 3, gv);



    /* merge gradient file - vs*/
    sprintf(modfile,"%s_gradient_vs_shot.%i.%i.%i",vinv->GRADIENT,ishot,gv->POS[1],gv->POS[2]);
	FP6=fopen(modfile,"wb");

	for (i=1;i<=gv->NX;i=i+gv->IDX){
		for (j=1;j<=gv->NY;j=j+gv->IDY){
			fwrite(&minv->gradVs_shot[j][i],sizeof(float),1,FP6);
		}
	}

	fclose(FP6);

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(modfile,"%s_gradient_vs_shot.%i",vinv->GRADIENT,ishot);
	if (gv->MPID==0) mergemod(modfile, 3, gv);

    writemod(modfile, minv->gradVs_shot, 3, gv);


    /* merge gradient filer - rho'*/
    sprintf(modfile,"%s_gradient_rho_shot.%i.%i.%i",vinv->GRADIENT,ishot,gv->POS[1],gv->POS[2]);
	FP6=fopen(modfile,"wb");

	for (i=1;i<=gv->NX;i=i+gv->IDX){
		for (j=1;j<=gv->NY;j=j+gv->IDY){
			fwrite(&minv->gradRho_shot[j][i],sizeof(float),1,FP6);
		}
	}

	fclose(FP6);

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(modfile,"%s_gradient_rho_shot.%i",vinv->GRADIENT,ishot);
	if (gv->MPID==0) mergemod(modfile, 3, gv);

    writemod(modfile, minv->gradRho_shot, 3, gv);

    /* output of gradient summed up for all shots */
    if (ishot == nshots){

        /* merge gradient file - vp*/
        sprintf(modfile,"%s_gradient.vp.%i.%i",vinv->GRADIENT,gv->POS[1],gv->POS[2]);
        FP6=fopen(modfile,"wb");

        for (i=1;i<=gv->NX;i=i+gv->IDX){
            for (j=1;j<=gv->NY;j=j+gv->IDY){
                fwrite(&minv->gradVp[j][i],sizeof(float),1,FP6);
            }
        }

        fclose(FP6);

        MPI_Barrier(MPI_COMM_WORLD);

        sprintf(modfile,"%s_gradient.vp",vinv->GRADIENT);
        if (gv->MPID==0) mergemod(modfile, 3, gv);

        writemod(modfile, minv->gradVp, 3, gv);



        /* merge gradient file - vs*/
        sprintf(modfile,"%s_gradient.vs.%i.%i",vinv->GRADIENT,gv->POS[1],gv->POS[2]);
        FP6=fopen(modfile,"wb");

        for (i=1;i<=gv->NX;i=i+gv->IDX){
            for (j=1;j<=gv->NY;j=j+gv->IDY){
                fwrite(&minv->gradVs[j][i],sizeof(float),1,FP6);
            }
        }

        fclose(FP6);

        MPI_Barrier(MPI_COMM_WORLD);

        sprintf(modfile,"%s_gradient.vs",vinv->GRADIENT);
        if (gv->MPID==0) mergemod(modfile, 3, gv);

        writemod(modfile, minv->gradVs, 3, gv);


        /* merge gradient filer - rho'*/
        sprintf(modfile,"%s_gradient.rho.%i.%i",vinv->GRADIENT,gv->POS[1],gv->POS[2]);
        FP6=fopen(modfile,"wb");

        for (i=1;i<=gv->NX;i=i+gv->IDX){
            for (j=1;j<=gv->NY;j=j+gv->IDY){
                fwrite(&minv->gradRho[j][i],sizeof(float),1,FP6);
            }
        }

        fclose(FP6);

        MPI_Barrier(MPI_COMM_WORLD);

        sprintf(modfile,"%s_gradient.rho",vinv->GRADIENT);
        if (gv->MPID==0) mergemod(modfile, 3, gv);

        writemod(modfile, minv->gradRho, 3, gv);

    }

}
