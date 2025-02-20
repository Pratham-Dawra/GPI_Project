
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* $Id: PML_pro.c 819 2015-04-17 11:07:06Z tmetz $ */

/*
* Define damping profiles for CPML boundary condition
* This C-PML implementation is adapted from the 2nd order isotropic CPML code by Dimitri Komatitsch and based in part on formulas given in Roden and Gedney (2000). 
* Additionally the code is based on the following references:
* 
* @ARTICLE{KoMa07,
* author = {Dimitri Komatitsch and Roland Martin},
* title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
*          at grazing incidence for the seismic wave equation},
* journal = {Geophysics},
* year = {2007},
* volume = {72},
* number = {5},
* pages = {SM155-SM167},
* doi = {10.1190/1.2757586}}
*
* @ARTICLE{MaKoEz08,
* author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
* title = {An unsplit convolutional perfectly matched layer improved at grazing
* incidence for seismic wave equation in poroelastic media},
* journal = {Geophysics},
* year = {2008},
* volume = {73},
* pages = {T51-T61},
* number = {4},
* doi = {10.1190/1.2939484}}
*
* @ARTICLE{MaKo09,
* author = {Roland Martin and Dimitri Komatitsch},
* title = {An unsplit convolutional perfectly matched layer technique improved
* at grazing incidence for the viscoelastic wave equation},
* journal = {Geophysical Journal International},
* year = {2009},
* volume = {179},
* pages = {333-344},
* number = {1},
* doi = {10.1111/j.1365-246X.2009.04278.x}}
*
* @ARTICLE{MaKoGe08,
* author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
* title = {A variational formulation of a stabilized unsplit convolutional perfectly
* matched layer for the isotropic or anisotropic seismic wave equation},
* journal = {Computer Modeling in Engineering and Sciences},
* year = {2008},
* volume = {37},
* pages = {274-304},
* number = {3}}
*
* The original CPML technique for Maxwell's equations is described in:
*
* @ARTICLE{RoGe00,
* author = {J. A. Roden and S. D. Gedney},
* title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
*          of the {CFS}-{PML} for Arbitrary Media},
* journal = {Microwave and Optical Technology Letters},
* year = {2000},
* volume = {27},
* number = {5},
* pages = {334-339},
* doi = {10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A}}
*/

#include "fd.h"

void PML_pro(MemModel * mpm, GlobVar * gv)
{
    int h;

    const float alpha_max_PML = 2.0 * PI * (gv->FPML / 2.0);    /* from festa and Vilotte */

    float thickness_PML_x, thickness_PML_y, xoriginleft, xoriginright, yoriginbottom, yorigintop;
    float Rcoef, d0_x, d0_y, xval, yval, abscissa_in_PML, abscissa_normalized;

    /* define profile of absorption in PML region */

    /* thickness of the PML layer in meters */
    thickness_PML_x = (float)gv->FW * gv->DH;
    thickness_PML_y = (float)gv->FW * gv->DH;

    /* reflection coefficient (INRIA report section 6.1) */
    Rcoef = 0.001;

    /* compute d0 from INRIA report section 6.1 */
    d0_x = -(gv->NPOWER + 1) * gv->VPPML * log(Rcoef) / (2.0 * thickness_PML_x);
    d0_y = -(gv->NPOWER + 1) * gv->VPPML * log(Rcoef) / (2.0 * thickness_PML_y);

    /* damping in the X direction */
    /* -------------------------- */

    /* origin of the PML layer (position of right edge minus thickness, in meters) */
    xoriginleft = thickness_PML_x;
    xoriginright = (gv->NXG - 1) * gv->DH - thickness_PML_x;

    /* left boundary */

    for (int i = 1; i <= gv->FW; i++) {

        mpm->K_x[i] = 1.0;
        mpm->K_x_half[i] = 1.0;
        xval = gv->DH * (i - 1);

        /* define damping profile at the grid points */
        abscissa_in_PML = xoriginleft - xval;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            mpm->d_x[i] = d0_x * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_x[i] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_x[i] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* define damping profile at half the grid points */
        abscissa_in_PML = xoriginleft - (xval + gv->DH / 2.0);

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            mpm->d_x_half[i] = d0_x * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_x_half[i] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_x_half[i] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* just in case, for -5 at the end */
        if (mpm->alpha_prime_x[i] < 0.0) {
            mpm->alpha_prime_x[i] = 0.0;
        }
        if (mpm->alpha_prime_x_half[i] < 0.0) {
            mpm->alpha_prime_x_half[i] = 0.0;
        }

        mpm->b_x[i] = exp(-(mpm->d_x[i] / mpm->K_x[i] + mpm->alpha_prime_x[i]) * gv->DT);
        mpm->b_x_half[i] = exp(-(mpm->d_x_half[i] / mpm->K_x_half[i] + mpm->alpha_prime_x_half[i]) * gv->DT);

        /* avoid division by zero outside the PML */
        if (fabsf(mpm->d_x[i]) > 1.0e-6) {
            mpm->a_x[i] =
                mpm->d_x[i] * (mpm->b_x[i] - 1.0) / (mpm->K_x[i] * (mpm->d_x[i] + mpm->K_x[i] * mpm->alpha_prime_x[i]));
        }
        if (fabsf(mpm->d_x_half[i]) > 1.0e-6) {
            mpm->a_x_half[i] =
                mpm->d_x_half[i] * (mpm->b_x_half[i] -
                                    1.0) / (mpm->K_x_half[i] * (mpm->d_x_half[i] +
                                                                mpm->K_x_half[i] * mpm->alpha_prime_x_half[i]));
        }

    }                           /* end of left boundary */

    /* right boundary */

    for (int i = gv->NXG - gv->FW + 1; i <= gv->NXG; i++) {

        h = i - gv->NXG + 2 * gv->FW;

        mpm->K_x[h] = 1.0;
        mpm->K_x_half[h] = 1.0;
        xval = gv->DH * (i - 1);

        /* define damping profile at the grid points */
        abscissa_in_PML = xval - xoriginright;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            mpm->d_x[h] = d0_x * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_x[h] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_x[h] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* define damping profile at half the grid points */
        abscissa_in_PML = xval + gv->DH / 2.0 - xoriginright;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            mpm->d_x_half[h] = d0_x * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_x_half[h] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_x_half[h] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* just in case, for -5 at the end */
        if (mpm->alpha_prime_x[h] < 0.0) {
            mpm->alpha_prime_x[h] = 0.0;
        }
        if (mpm->alpha_prime_x_half[h] < 0.0) {
            mpm->alpha_prime_x_half[h] = 0.0;
        }

        mpm->b_x[h] = exp(-(mpm->d_x[h] / mpm->K_x[h] + mpm->alpha_prime_x[h]) * gv->DT);
        mpm->b_x_half[h] = exp(-(mpm->d_x_half[h] / mpm->K_x_half[h] + mpm->alpha_prime_x_half[h]) * gv->DT);

        /* avoid division by zero outside the PML */
        if (fabsf(mpm->d_x[h]) > 1.0e-6) {
            mpm->a_x[h] =
                mpm->d_x[h] * (mpm->b_x[h] - 1.0) / (mpm->K_x[h] * (mpm->d_x[h] + mpm->K_x[h] * mpm->alpha_prime_x[h]));
        }
        if (fabsf(mpm->d_x_half[h]) > 1.0e-6) {
            mpm->a_x_half[h] =
                mpm->d_x_half[h] * (mpm->b_x_half[h] -
                                    1.0) / (mpm->K_x_half[h] * (mpm->d_x_half[h] +
                                                                mpm->K_x_half[h] * mpm->alpha_prime_x_half[h]));
        }

    }                           /* end of right boundary */

    /* damping in the Y direction */
    /* -------------------------- */

    /* origin of the PML layer (position of right edge minus thickness, in meters) */
    yoriginbottom = thickness_PML_y;
    yorigintop = (gv->NYG - 1) * gv->DH - thickness_PML_y;

    for (int i = 1; i <= gv->FW; i++) {

        mpm->K_y[i] = 1.0;
        mpm->K_y_half[i] = 1.0;
        yval = gv->DH * (i - 1);

        /* left boundary */

        /* define damping profile at the grid points */
        abscissa_in_PML = yoriginbottom - yval;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_y;
            mpm->d_y[i] = d0_y * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_y[i] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_y[i] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* define damping profile at half the grid points */
        abscissa_in_PML = yoriginbottom - (yval + gv->DH / 2.0);

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_y;
            mpm->d_y_half[i] = d0_y * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_y_half[i] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_y_half[i] = alpha_max_PML * (1.0 - abscissa_normalized);
        }
        mpm->b_y[i] = exp(-(mpm->d_y[i] / mpm->K_y[i] + mpm->alpha_prime_y[i]) * gv->DT);
        mpm->b_y_half[i] = exp(-(mpm->d_y_half[i] / mpm->K_y_half[i] + mpm->alpha_prime_y_half[i]) * gv->DT);

        /* avoid division by zero outside the PML */
        if (fabsf(mpm->d_y[i]) > 1.0e-6) {
            mpm->a_y[i] =
                mpm->d_y[i] * (mpm->b_y[i] - 1.0) / (mpm->K_y[i] * (mpm->d_y[i] + mpm->K_y[i] * mpm->alpha_prime_y[i]));
        }
        if (fabsf(mpm->d_y_half[i]) > 1.0e-6) {
            mpm->a_y_half[i] =
                mpm->d_y_half[i] * (mpm->b_y_half[i] -
                                    1.0) / (mpm->K_y_half[i] * (mpm->d_y_half[i] +
                                                                mpm->K_y_half[i] * mpm->alpha_prime_y_half[i]));
        }

    }                           /* end of left boundary */

    /* top boundary */
    for (int i = gv->NYG - gv->FW + 1; i <= gv->NYG; i++) {

        h = i - gv->NYG + 2 * gv->FW;

        mpm->K_y[h] = 1.0;
        mpm->K_y_half[h] = 1.0;
        yval = gv->DH * (i - 1);

        /* define damping profile at the grid points */
        abscissa_in_PML = yval - yorigintop;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_y;
            mpm->d_y[h] = d0_y * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_y[h] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_y[h] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        /* define damping profile at half the grid points */
        abscissa_in_PML = yval + gv->DH / 2.0 - yorigintop;

        if (abscissa_in_PML >= 0.0) {
            abscissa_normalized = abscissa_in_PML / thickness_PML_y;
            mpm->d_y_half[h] = d0_y * pow(abscissa_normalized, gv->NPOWER);

            /* this taken from Gedney page 8.2 */
            mpm->K_y_half[h] = 1.0 + (gv->K_MAX_CPML - 1.0) * pow(abscissa_normalized, gv->NPOWER);
            mpm->alpha_prime_y_half[h] = alpha_max_PML * (1.0 - abscissa_normalized);
        }

        mpm->b_y[h] = exp(-(mpm->d_y[h] / mpm->K_y[h] + mpm->alpha_prime_y[h]) * gv->DT);
        mpm->b_y_half[h] = exp(-(mpm->d_y_half[h] / mpm->K_y_half[h] + mpm->alpha_prime_y_half[h]) * gv->DT);

        /* avoid division by zero outside the PML */
        if (fabsf(mpm->d_y[h]) > 1.0e-6) {
            mpm->a_y[h] =
                mpm->d_y[h] * (mpm->b_y[h] - 1.0) / (mpm->K_y[h] * (mpm->d_y[h] + mpm->K_y[h] * mpm->alpha_prime_y[h]));
        }
        if (fabsf(mpm->d_y_half[h]) > 1.0e-6) {
            mpm->a_y_half[h] =
                mpm->d_y_half[h] * (mpm->b_y_half[h] -
                                    1.0) / (mpm->K_y_half[h] * (mpm->d_y_half[h] +
                                                                mpm->K_y_half[h] * mpm->alpha_prime_y_half[h]));
        }
    }                           /* end of top boundary */
}
