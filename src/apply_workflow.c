
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

/*------------------------------------------------------------------------
 *  Apply workflow
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "macros.h"
#include "logging.h"

void apply_workflow(int iter, GlobVar *gv, GlobVarInv *vinv)
{

    /******************/
    /* Apply Workflow */
    /******************/

    /* Inversion of material parameter */
    if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][1] != -1) {
        if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][1] == 1) {
            if (vinv->INV_VS_ITER > iter)
                vinv->INV_VS_ITER = iter;
        } else {
            /* detect change and reset LBFGS */
            if (vinv->INV_VS_ITER < iter)
                vinv->LBFGS_ITER_START = iter;
            vinv->INV_VS_ITER = iter + 10;
        }
    }

    if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][2] != -1) {
        if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][2] == 1) {
            if (vinv->INV_VP_ITER > iter)
                vinv->INV_VP_ITER = iter;
        } else {
            /* detect change and reset LBFGS */
            if (vinv->INV_VP_ITER < iter)
                vinv->LBFGS_ITER_START = iter;
            vinv->INV_VP_ITER = iter + 10;
        }
    }

    if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][3] != -1) {
        if (vinv->WORKFLOW[vinv->WORKFLOW_STAGE][3] == 1) {
            if (vinv->INV_RHO_ITER > iter)
                vinv->INV_RHO_ITER = iter;
        } else {
            /* detect change and reset LBFGS */
            if (vinv->INV_RHO_ITER < iter)
                vinv->LBFGS_ITER_START = iter;
            vinv->INV_RHO_ITER = iter + 10;
        }
    }

    /* Abort criterium */
    vinv->PRO = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][4];

    /* Frequency filtering  */
    vinv->TIME_FILT = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][5];

    if (vinv->TIME_FILT > 0) {
        if (vinv->F_HIGH_PASS != vinv->WORKFLOW[vinv->WORKFLOW_STAGE][6])
            vinv->LBFGS_ITER_START = iter;
        vinv->F_HIGH_PASS = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][6];

        if (*(vinv->F_LOW_PASS) != vinv->WORKFLOW[vinv->WORKFLOW_STAGE][7])
            vinv->LBFGS_ITER_START = iter;
        *(vinv->F_LOW_PASS) = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][7];
    }

    /* Approx. Hessian  */
    /*if (vinv->EPRECOND == 0 && vinv->WORKFLOW[vinv->WORKFLOW_STAGE][8] != 0) {
     * if (gv->MPID == 0)
     * log_warn("WARNING: EPRECOND have to be set >0 in JSON (if so, ignore this message).\n");
     * } */

    vinv->EPRECOND = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][8];
    vinv->EPSILON_WE = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][9];
    vinv->TIMEWIN = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][10];
    vinv->GAMMA = vinv->WORKFLOW[vinv->WORKFLOW_STAGE][11];

    if (vinv->LBFGS_ITER_START == iter && vinv->GRAD_METHOD == 2) {
        if (gv->MPID == 0)
            log_info("L-BFGS will be used from iteration %d on.\n", vinv->LBFGS_ITER_START + 1);
    }

    /* Print current workflow */
    if (gv->MPID == 0) {
        log_info("---------- Applying Workflow -----------\n");
        log_info("INV_VS_ITER: %d\n", vinv->INV_VS_ITER);
        log_info("INV_VP_ITER: %d\n", vinv->INV_VP_ITER);
        log_info("INV_RHO_ITER: %d\n", vinv->INV_RHO_ITER);
        log_info("PRO: %f\n", vinv->PRO);
        log_info("TIME_FILT: %d\n", vinv->TIME_FILT);
        log_info("F_HIGH_PASS: %f\n", vinv->F_HIGH_PASS);
        log_info("F_LOW_PASS: %f\n", *(vinv->F_LOW_PASS));
        log_info("EPRECOND: %d\n", vinv->EPRECOND);
        log_info("EPSILON_WE: %f\n", vinv->EPSILON_WE);
        log_info("TIMEWIN: %d\n", vinv->TIMEWIN);
        log_info("GAMMA: %f\n", vinv->GAMMA);
    }
}
