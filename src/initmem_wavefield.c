
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

/* -------------------------------------------------------------
 * This is function initmem_elastic.
 * Initialising memory for wavefield variables.
 * -------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void initmem_wavefield(MemWavefield *mpw, GlobVar *gv)
{
    switch (gv->WEQ) {
      case AC_ISO:             /* acoustic */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = mpw->psxx;
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->dummy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          break;
      case EL_ISO:             /* elastic */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          break;
      case VEL_ISO:            /* viscoelastic */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pr = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pp = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pq = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          break;
      case EL_VTI:             /* elastic VTI */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          break;
      case EL_TTI:             /* elastic TTI */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          break;
      case VEL_VTI:            /* elastic VTI */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pr = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pp = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pq = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          break;
      case VEL_TTI:            /* elastic VTI */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->psyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pr = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pp = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          mpw->pq = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          break;
      case VAC_ISO:            /* viscoacoustic */
          mpw->psxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
          mpw->pr = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
          break;
      default:
          log_fatal("Unknown WEQ.\n");
    }

    mpw->pvxx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    mpw->pvyy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    mpw->pvyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    mpw->pvxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    if (gv->FDORDER_TIME == 4) {
        mpw->vxx_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxx_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxx_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxx_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpw->vyy_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyy_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyy_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyy_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpw->vxy_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxy_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxy_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vxy_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpw->vyx_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyx_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyx_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->vyx_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpw->svx_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svx_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svx_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svx_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpw->svy_1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svy_2 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svy_3 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpw->svy_4 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->L > 0 && gv->FDORDER_TIME == 4) {
        mpw->pr_2 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pr_3 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pr_4 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        mpw->pp_2 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pp_3 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pp_4 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        mpw->pq_2 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pq_3 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpw->pq_4 = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
    }

    if (gv->ABS_TYPE == 1) {    /* PML */
        mpw->psi_sxx_x = matrix(1, gv->NY, 1, 2 * gv->FW);
        mpw->psi_syy_y = matrix(1, 2 * gv->FW, 1, gv->NX);
        mpw->psi_sxy_y = matrix(1, 2 * gv->FW, 1, gv->NX);
        mpw->psi_sxy_x = matrix(1, gv->NY, 1, 2 * gv->FW);

        mpw->psi_vxx = matrix(1, gv->NY, 1, 2 * gv->FW);
        mpw->psi_vyy = matrix(1, 2 * gv->FW, 1, gv->NX);
        mpw->psi_vxy = matrix(1, 2 * gv->FW, 1, gv->NX);
        mpw->psi_vyx = matrix(1, gv->NY, 1, 2 * gv->FW);

        mpw->psi_vxxs = matrix(1, gv->NY, 1, 2 * gv->FW);   /* For surface_elastic(visc).c */
    }

    /* memory allocation for buffer arrays in which the wavefield
     * information to be exchanged between neighboring PEs is stored */
    mpw->bufferlef_to_rig = matrix(1, gv->NY, 1, gv->FDORDER);
    mpw->bufferrig_to_lef = matrix(1, gv->NY, 1, gv->FDORDER);
    mpw->buffertop_to_bot = matrix(1, gv->NX, 1, gv->FDORDER);
    mpw->bufferbot_to_top = matrix(1, gv->NX, 1, gv->FDORDER);
}
