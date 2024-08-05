
/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------
 * CPML update functions (ABS=1) for the stress and particle velocitiy components used
 * in the update functions of the CPML-boundary area.
 * -------------------------------------------------------------------------------*/

#include "fd.h"

/* CPML Functions for update_s ---------------------------------------------------*/

void cpml_update_s_x(int i, int j, int h1, int h2, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_vxx[h2][h1] = mpm->b_x[h1] * mpw->psi_vxx[h2][h1] + mpm->a_x[h1] * (mpw->pvxx[j][i]);
    mpw->pvxx[j][i] = (mpw->pvxx[j][i]) / mpm->K_x[h1] + mpw->psi_vxx[h2][h1];

    mpw->psi_vyx[h2][h1] = mpm->b_x_half[h1] * mpw->psi_vyx[h2][h1] + mpm->a_x_half[h1] * (mpw->pvyx[j][i]);
    mpw->pvyx[j][i] = (mpw->pvyx[j][i]) / mpm->K_x_half[h1] + mpw->psi_vyx[h2][h1];
}

void cpml_update_s_y(int i, int j, int h1, int h2, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_vyy[h2][h1] = mpm->b_y[h2] * mpw->psi_vyy[h2][h1] + mpm->a_y[h2] * (mpw->pvyy[j][i]);
    mpw->pvyy[j][i] = (mpw->pvyy[j][i]) / mpm->K_y[h2] + mpw->psi_vyy[h2][h1];

    mpw->psi_vxy[h2][h1] = mpm->b_y_half[h2] * mpw->psi_vxy[h2][h1] + mpm->a_y_half[h2] * (mpw->pvxy[j][i]);
    mpw->pvxy[j][i] = (mpw->pvxy[j][i]) / mpm->K_y_half[h2] + mpw->psi_vxy[h2][h1];
}

/* CPML Functions for update_v ---------------------------------------------------*/

void cpml_update_v_x(int i, int j, float *sxx_x, float *sxy_x, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_sxx_x[j][i] = mpm->b_x_half[i] * mpw->psi_sxx_x[j][i] + mpm->a_x_half[i] * (*sxx_x);
    *sxx_x = (*sxx_x) / mpm->K_x_half[i] + mpw->psi_sxx_x[j][i];

    mpw->psi_sxy_x[j][i] = mpm->b_x[i] * mpw->psi_sxy_x[j][i] + mpm->a_x[i] * (*sxy_x);
    *sxy_x = (*sxy_x) / mpm->K_x[i] + mpw->psi_sxy_x[j][i];
}

void cpml_update_v_y(int i, int j, float *sxy_y, float *syy_y, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_syy_y[j][i] = mpm->b_y_half[j] * mpw->psi_syy_y[j][i] + mpm->a_y_half[j] * (*syy_y);
    *syy_y = (*syy_y) / mpm->K_y_half[j] + mpw->psi_syy_y[j][i];

    mpw->psi_sxy_y[j][i] = mpm->b_y[j] * mpw->psi_sxy_y[j][i] + mpm->a_y[j] * (*sxy_y);
    *sxy_y = (*sxy_y) / mpm->K_y[j] + mpw->psi_sxy_y[j][i];
}

/* acoustic cases ---------------------------------------------------*/

void cpml_update_s_x_ac(int i, int j, int h1, int h2, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_vxx[h2][h1] = mpm->b_x[h1] * mpw->psi_vxx[h2][h1] + mpm->a_x[h1] * (mpw->pvxx[j][i]);
    mpw->pvxx[j][i] = (mpw->pvxx[j][i]) / mpm->K_x[h1] + mpw->psi_vxx[h2][h1];
}

void cpml_update_s_y_ac(int i, int j, int h1, int h2, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_vyy[h2][h1] = mpm->b_y[h2] * mpw->psi_vyy[h2][h1] + mpm->a_y[h2] * (mpw->pvyy[j][i]);
    mpw->pvyy[j][i] = (mpw->pvyy[j][i]) / mpm->K_y[h2] + mpw->psi_vyy[h2][h1];
}

/* CPML Functions for update_v ---------------------------------------------------*/

void cpml_update_v_x_ac(int i, int j, float *sxx_x, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_sxx_x[j][i] = mpm->b_x_half[i] * mpw->psi_sxx_x[j][i] + mpm->a_x_half[i] * (*sxx_x);
    *sxx_x = (*sxx_x) / mpm->K_x_half[i] + mpw->psi_sxx_x[j][i];

}

void cpml_update_v_y_ac(int i, int j, float *syy_y, MemModel * mpm, MemWavefield * mpw)
{
    mpw->psi_syy_y[j][i] = mpm->b_y_half[j] * mpw->psi_syy_y[j][i] + mpm->a_y_half[j] * (*syy_y);
    *syy_y = (*syy_y) / mpm->K_y_half[j] + mpw->psi_syy_y[j][i];

}
