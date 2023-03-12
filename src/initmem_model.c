
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
 * This is function initmem_model.
 * Initialising memory for model variables.
 * -----------------------------------------------------------*/

#include "fd.h"

void initmem_model(MemModel *mpm, GlobVar *gv)
{

    /* static (model) arrays (isotropic elastic + viscoelastic) */
    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {
        mpm->prho = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->prip = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->prjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ppi = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->puipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->absorb_coeff = matrix(1, gv->NY, 1, gv->NX);
    }

    /* static (model) arrays (viscoelastic) */
    if (gv->WEQ == VEL_ISO) {   /*viscoelastic  isotropic wave equation */
        mpm->dip = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->e = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->ptaus = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptausipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptaup = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->fipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->f = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->g = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->peta = vector(1, gv->L);
        mpm->bip = vector(1, gv->L);
        mpm->bjm = vector(1, gv->L);
        mpm->cip = vector(1, gv->L);
        mpm->cjm = vector(1, gv->L);
    }

    if (gv->WEQ == EL_VTI) {    /*elastic VTI wave equation */
        mpm->pc11 = mpm->ppi;
        mpm->pc33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc55 = mpm->pu;
        mpm->pc55ipjp = mpm->puipjp;
    }

    if (gv->WEQ == VEL_VTI) {   /*viscoelastic VTI wave equation */
        mpm->pc11 = mpm->ppi;
        mpm->pc33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc55 = mpm->pu;
        mpm->pc55ipjp = mpm->puipjp;

        mpm->ptau11 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau15 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau55 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau55ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->pc55ipjpd = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc13d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc33d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc11d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        mpm->pc55ipjpu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc11u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc33u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->peta = vector(1, gv->L);
        mpm->bip = vector(1, gv->L);
        mpm->cip = vector(1, gv->L);
    }

    if (gv->WEQ == EL_TTI) {    /*elastic TTI wave equation */
        mpm->pc11 = mpm->ppi;
        mpm->pc33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc55 = mpm->pu;
        mpm->pc55ipjp = mpm->puipjp;
        mpm->pc15 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc15ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    if (gv->WEQ == VEL_TTI) {   /*viscoelastic TTI wave equation */
        mpm->pc11 = mpm->ppi;
        mpm->pc33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc55 = mpm->pu;
        mpm->pc15 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->pc55ipjp = mpm->puipjp;
        mpm->pc15ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->ptau11 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau33 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau13 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau55 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau15 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau35 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau55ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau15ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->ptau35ipjp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->pc11d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc33d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc13d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc55d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc15d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc35d = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc55ipjpd = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc15ipjpd = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);
        mpm->pc35ipjpd = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, gv->L);

        mpm->pc11u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc33u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc13u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc55u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc15u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->pc55ipjpu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc15ipjpu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        mpm->pc35ipjpu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        mpm->peta = vector(1, gv->L);
        mpm->bip = vector(1, gv->L);
        mpm->cip = vector(1, gv->L);
    }

    if (gv->ABS_TYPE == 1) {    /* PML */
        mpm->d_x = vector(1, 2 * gv->FW);
        mpm->K_x = vector(1, 2 * gv->FW);
        mpm->alpha_prime_x = vector(1, 2 * gv->FW);
        mpm->a_x = vector(1, 2 * gv->FW);
        mpm->b_x = vector(1, 2 * gv->FW);

        mpm->d_x_half = vector(1, 2 * gv->FW);
        mpm->K_x_half = vector(1, 2 * gv->FW);
        mpm->alpha_prime_x_half = vector(1, 2 * gv->FW);
        mpm->a_x_half = vector(1, 2 * gv->FW);
        mpm->b_x_half = vector(1, 2 * gv->FW);

        mpm->d_y = vector(1, 2 * gv->FW);
        mpm->K_y = vector(1, 2 * gv->FW);
        mpm->alpha_prime_y = vector(1, 2 * gv->FW);
        mpm->a_y = vector(1, 2 * gv->FW);
        mpm->b_y = vector(1, 2 * gv->FW);

        mpm->d_y_half = vector(1, 2 * gv->FW);
        mpm->K_y_half = vector(1, 2 * gv->FW);
        mpm->alpha_prime_y_half = vector(1, 2 * gv->FW);
        mpm->a_y_half = vector(1, 2 * gv->FW);
        mpm->b_y_half = vector(1, 2 * gv->FW);
    }
}
