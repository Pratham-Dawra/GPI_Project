

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
 * Calculate Convolution of Forward and Backward model for every DTINV time step
 *-------------------------------------------------------------------------------*/
#include "fd.h"
#include "logging.h"

//float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int sws,
//int swstestshot, int ishot, int iter, AcqVar *acq, GlobVarInv *vinv

void calc_conv(int hin, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    int IDXI = 1, IDYI = 1;
    float muss, lamss;
    int i, j;

    /*-------------------------------------------------*/
    /* Calculate convolution  */

    /*-------------------------------------------------*/

    for (j = 1; j <= gv->NY; j = j + IDYI) {
        for (i = 1; i <= gv->NX; i = i + IDXI) {

            /* get gradient for lambda/mu/rho parameterization */

            muss = minv->Vs0[j][i] * minv->Vs0[j][i] * minv->Rho0[j][i];
            lamss = minv->Vp0[j][i] * minv->Vp0[j][i] * minv->Rho0[j][i] - 2.0 * muss;

            minv->gradRhos_shot[j][i] += - gv->DT *
                ((minv->pvxp1[j][i] * minv->forward_prop_rho_x[j][i][vinv->NTDTINV - hin + 1]) +
                (minv->pvyp1[j][i] * minv->forward_prop_rho_y[j][i][vinv->NTDTINV - hin + 1]));

            minv->gradLam_shot[j][i] += - gv->DT * (1.0/(4.0 * (lamss+muss) * (lamss+muss))) *
                ((minv->forward_prop_x[j][i][vinv->NTDTINV - hin + 1] +
                 minv->forward_prop_y[j][i][vinv->NTDTINV - hin + 1]) * (mpw->psxx[j][i] + mpw->psyy[j][i]));

            if (mpm->pu[j][i] > 0.0) {
                minv->gradMu_shot[j][i] += - gv->DT *
                    (((1.0 / (muss * muss)) * (minv->forward_prop_u[j][i][vinv->NTDTINV - hin + 1] *
                    mpw->psxy[j][i])) +
                    ((1.0 / 4.0) *
                    ((minv->forward_prop_x[j][i][vinv->NTDTINV - hin + 1] +
                    minv->forward_prop_y[j][i][vinv->NTDTINV - hin + 1]) * (mpw->psxx[j][i] +
                    mpw->psyy[j][i])) / ((lamss + muss) * (lamss + muss))) +
                    ((1.0 / 4.0) *
                    ((minv->forward_prop_x[j][i][vinv->NTDTINV - hin + 1] -
                    minv->forward_prop_y[j][i][vinv->NTDTINV - hin + 1]) *
                    (mpw->psxx[j][i] - mpw->psyy[j][i])) / (muss * muss)));
            }

        }
    }
}
