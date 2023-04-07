
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

void initmem_fwi(MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    /* Wavefield allocations */
    minv->ux = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->uy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->uxy = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->uyx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->uttx = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->utty = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->pvxp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->pvyp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->pvxm1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->pvym1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    /* Model allocations */
    minv->Vp0 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->Vs0 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->Rho0 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->prhonp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->pripnp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->prjpnp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->ppinp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->punp1 = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->vpmat = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    /* other allocations */
    if ((vinv->EPRECOND == 1) || (vinv->EPRECOND == 3)) {
        minv->We_sum = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        minv->Ws = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);  /* total energy of the source wavefield */
        minv->Wr = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);  /* total energy of the receiver wavefield */
        minv->We = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);  /* total energy of source agv->ND receiver wavefield */
    }

    if (vinv->GRAD_METHOD == 2) {
        /* Allocate memory for L-BFGS */
        minv->y_LBFGS = fmatrix(1, vinv->N_LBFGS, 1, LBFGS_NPAR * gv->NX * gv->NY);
        minv->s_LBFGS = fmatrix(1, vinv->N_LBFGS, 1, LBFGS_NPAR * gv->NX * gv->NY);
        minv->rho_LBFGS = vector(1, vinv->N_LBFGS);
    }

    minv->waveconv = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_lam = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_shot = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->waveconvtmp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->wcpart = matrix(1, 3, 1, 3);
    minv->wavejac = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->forward_prop_x = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, vinv->NTDTINV);
    minv->forward_prop_y = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, vinv->NTDTINV);

    minv->gradg = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->gradp = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->forward_prop_rho_x = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, vinv->NTDTINV);
    minv->forward_prop_rho_y = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, vinv->NTDTINV);

    minv->gradg_rho = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->gradp_rho = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_rho = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_rho_s = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_rho_shot = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    if (vinv->WOLFE_CONDITION) {

        minv->waveconv_old = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        minv->waveconv_u_old = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        minv->waveconv_rho_old = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

        minv->waveconv_up = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        minv->waveconv_u_up = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
        minv->waveconv_rho_up = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    }

    minv->forward_prop_u = f3tensor(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND, 1, vinv->NTDTINV);
    minv->gradg_u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->gradp_u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_u = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_mu = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);
    minv->waveconv_u_shot = matrix(-gv->ND + 1, gv->NY + gv->ND, -gv->ND + 1, gv->NX + gv->ND);

    minv->taper_coeff = matrix(1, gv->NY, 1, gv->NX);

    /* Allocate memory to save full seismograms */
    switch (gv->SEISMO) {
      case 1:                  /* particle velocities only */
          minv->fulldata_vx = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_vy = matrix(1, gv->NTRG, 1, gv->NT);
          break;
      case 2:                  /* pressure only */
          minv->fulldata_p = matrix(1, gv->NTRG, 1, gv->NT);
          break;
      case 3:                  /* curl and div only */
          minv->fulldata_div = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_curl = matrix(1, gv->NTRG, 1, gv->NT);
          break;
      case 4:                  /* everything */
          minv->fulldata_vx = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_vy = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_p = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_div = matrix(1, gv->NTRG, 1, gv->NT);
          minv->fulldata_curl = matrix(1, gv->NTRG, 1, gv->NT);
          break;
    }

    minv->sectionpnp1 = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionpn = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionpdata = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionpdiff = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionpdiffold = matrix(1, gv->NTR, 1, gv->NS);

    minv->sectionvxdata = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionvxdiff = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionvxdiffold = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionvydata = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionvydiff = matrix(1, gv->NTR, 1, gv->NS);
    minv->sectionvydiffold = matrix(1, gv->NTR, 1, gv->NS);

    /* Memory for seismic data */
    minv->sectionread = matrix(1, gv->NTRG, 1, gv->NS);

    /* Memory for inversion for source time function */
    if ((vinv->INV_STF == 1) || (vinv->TIME_FILT == 1) || (vinv->TIME_FILT == 2)) {
        minv->source_time_function = vector(1, gv->NT);
        minv->sectionvy_conv = matrix(1, gv->NTRG, 1, gv->NT);
        minv->sectionvy_obs = matrix(1, gv->NTRG, 1, gv->NT);
        minv->sectionvx_conv = matrix(1, gv->NTRG, 1, gv->NT);
        minv->sectionvx_obs = matrix(1, gv->NTRG, 1, gv->NT);
        minv->sectionp_conv = matrix(1, gv->NTRG, 1, gv->NT);
        minv->sectionp_obs = matrix(1, gv->NTRG, 1, gv->NT);
    }

    /* memory of L2 norm */
    minv->L2t = vector(1, 4);
    minv->epst1 = vector(1, 3);
    minv->epst2 = vector(1, 3);
    minv->epst3 = vector(1, 3);
    minv->picked_times = vector(1, gv->NTR);

}
