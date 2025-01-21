
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

/* -----------------------------------------------------------
 * This is function freemem_model.
 * De-allocating memory of model variables.
 * -----------------------------------------------------------*/

#include "fd.h"

void freemem_model(MemModel *mpm, GlobVar *gv)
{

    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {
        free_matrix(mpm->prho, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->prip, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->prjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ppi, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pu, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->puipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->absorb_coeff, 1, gv->NY, 1, gv->NX);
    }

    if (gv->WEQ == VEL_ISO) {   /*viscoelastic wave equation */
        free_f3tensor(mpm->dip, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->e, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_matrix(mpm->ptaus, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptausipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptaup, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->fipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->f, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->g, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_vector(mpm->peta, 1, gv->L);
        free_vector(mpm->bip, 1, gv->L);
        free_vector(mpm->bjm, 1, gv->L);
        free_vector(mpm->cip, 1, gv->L);
        free_vector(mpm->cjm, 1, gv->L);
    }

    if (gv->WEQ == EL_VTI) {    /*elastic VTI wave equation */
        free_matrix(mpm->pc33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->WEQ == VEL_VTI) {   /*viscoelastic VTI wave equation */
        free_matrix(mpm->pc33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpm->ptau11, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau15, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau55, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau55ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_f3tensor(mpm->pc55ipjpd, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc13d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc33d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc11d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        free_matrix(mpm->pc55ipjpu, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc11u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc33u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_vector(mpm->peta, 1, gv->L);
        free_vector(mpm->bip, 1, gv->L);
        free_vector(mpm->cip, 1, gv->L);
    }

    if (gv->WEQ == EL_TTI) {    /*elastic TTI wave equation */
        free_matrix(mpm->pc33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc15, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc15ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->WEQ == VEL_TTI) {   /*viscoelastic TTI wave equation */
        free_matrix(mpm->pc33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc15, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpm->pc15ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpm->ptau11, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau33, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau13, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau55, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau15, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau35, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau55ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau15ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->ptau35ipjp, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_f3tensor(mpm->pc11d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc33d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc13d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc55d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc15d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc35d, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc55ipjpd, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc15ipjpd, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        free_f3tensor(mpm->pc35ipjpd, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        free_matrix(mpm->pc11u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc33u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc13u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc55u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc15u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35u, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_matrix(mpm->pc55ipjpu, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc15ipjpu, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        free_matrix(mpm->pc35ipjpu, -gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        free_vector(mpm->peta, 1, gv->L);
        free_vector(mpm->bip, 1, gv->L);
        free_vector(mpm->cip, 1, gv->L);
    }

    if (gv->ABS_TYPE == 1) {
        free_vector(mpm->d_x, 1, 2 * gv->FW);
        free_vector(mpm->K_x, 1, 2 * gv->FW);
        free_vector(mpm->alpha_prime_x, 1, 2 * gv->FW);
        free_vector(mpm->a_x, 1, 2 * gv->FW);
        free_vector(mpm->b_x, 1, 2 * gv->FW);

        free_vector(mpm->d_x_half, 1, 2 * gv->FW);
        free_vector(mpm->K_x_half, 1, 2 * gv->FW);
        free_vector(mpm->alpha_prime_x_half, 1, 2 * gv->FW);
        free_vector(mpm->a_x_half, 1, 2 * gv->FW);
        free_vector(mpm->b_x_half, 1, 2 * gv->FW);

        free_vector(mpm->d_y, 1, 2 * gv->FW);
        free_vector(mpm->K_y, 1, 2 * gv->FW);
        free_vector(mpm->alpha_prime_y, 1, 2 * gv->FW);
        free_vector(mpm->a_y, 1, 2 * gv->FW);
        free_vector(mpm->b_y, 1, 2 * gv->FW);

        free_vector(mpm->d_y_half, 1, 2 * gv->FW);
        free_vector(mpm->K_y_half, 1, 2 * gv->FW);
        free_vector(mpm->alpha_prime_y_half, 1, 2 * gv->FW);
        free_vector(mpm->a_y_half, 1, 2 * gv->FW);
        free_vector(mpm->b_y_half, 1, 2 * gv->FW);
    }
}
