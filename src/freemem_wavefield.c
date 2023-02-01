
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
 * This is function freemem_wavefield.
 * De-allocating memory of wavefield variables.
 *
 * -------------------------------------------------------------*/

#include "fd.h"

void freemem_wavefield(MemWavefield * mpw, GlobVar * gv)
{

    /* deallocation of memory */
    free_matrix(mpw->psxx, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    free_matrix(mpw->psxy, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    free_matrix(mpw->psyy, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    free_matrix(mpw->pvx, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    free_matrix(mpw->pvy, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    if (gv->FDORDER_TIME == 4) {
        free_matrix(mpw->vxx_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxx_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxx_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxx_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpw->vyy_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyy_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyy_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyy_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpw->vxy_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxy_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxy_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vxy_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpw->vyx_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyx_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyx_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->vyx_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpw->svx_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svx_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svx_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svx_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpw->svy_1, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svy_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svy_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->svy_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->L) {
        free_f3tensor(mpw->pr, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pq, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
    }

    if (gv->L > 0 && gv->FDORDER_TIME == 4) {
        free_f3tensor(mpw->pr_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pr_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pr_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        free_f3tensor(mpw->pp_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pp_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pp_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        free_f3tensor(mpw->pq_2, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pq_3, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpw->pq_4, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
    }

    /*elastic TTI wave equation */
    if (gv->WEQ == EL_TTI) {    /*elastic TTI wave equation */
        free_matrix(mpw->pvxx, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->pvyy, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->pvxy, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpw->pvyx, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->ABS_TYPE == 1) {
        free_matrix(mpw->psi_sxx_x, 1, gv->NY, 1, 2 * gv->FW);
        free_matrix(mpw->psi_syy_y, 1, 2 * gv->FW, 1, gv->NX);
        free_matrix(mpw->psi_sxy_x, 1, gv->NY, 1, 2 * gv->FW);
        free_matrix(mpw->psi_sxy_y, 1, 2 * gv->FW, 1, gv->NX);

        free_matrix(mpw->psi_vxx, 1, gv->NY, 1, 2 * gv->FW);
        free_matrix(mpw->psi_vyy, 1, 2 * gv->FW, 1, gv->NX);
        free_matrix(mpw->psi_vxy, 1, 2 * gv->FW, 1, gv->NX);
        free_matrix(mpw->psi_vyx, 1, gv->NY, 1, 2 * gv->FW);

        free_matrix(mpw->psi_vxxs, 1, gv->NY, 1, 2 * gv->FW);   /* For surface_elastic(visc).c */
    }
}
