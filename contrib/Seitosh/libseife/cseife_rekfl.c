/*! \file cseife_rekfl.c
 * \brief apply IIR filter coefficients (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * Copyright 1984 by Erhard Wielandt
 * This algorithm was part of seife.f. A current version of seife.f can be
 * obtained from http://www.software-for-seismometry.de/
 * 
 * ----
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * ----
 *
 * \date 28/06/2005
 * 
 * apply IIR filter coefficients (implementation)
 * 
 * REVISIONS and CHANGES 
 *  - 28/06/2005   V1.0   Thomas Forbriger
 *  - 11/07/2005   V1.1   support debug mode
 * 
 * ============================================================================
 */
#define TF_CSEIFE_REKFL_C_VERSION \
  "TF_CSEIFE_REKFL_C   V1.1"

#include <stdio.h>
#include "cseife.h"

/* subs/seife_rekfl.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

   the code was derived through f2c, but modified thereafter
*/

void seife_rekfl(double* x, double* y, int n, struct seife_IIRcoef coef)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int j;
    double xa, ya, xn, xaa, yaa;

/* report coefficients is requested */
    if (seife_global_flags.debug == 1)
    {
      fprintf(stderr, "DEBUG (seife_rekfl):"
              "%s=%-10.3g %s=%-10.3g %s=%-10.3g %s=%-10.3g %s=%-10.3g\n",
              "f0",coef.f0, "f1",coef.f1, "f2",coef.f2, "g1",coef.g1,
              "g2",coef.g2);
    }

/*  perform recursive filtering */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    xa = 0.;
    xaa = 0.;
    ya = 0.;
    yaa = 0.;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	xn = x[j];
	y[j] = coef.f0 * xn + coef.f1 * xa + coef.f2 * xaa + coef.g1 * ya +
        coef.g2 * yaa;
	xaa = xa;
	xa = xn;
	yaa = ya;
/* L1: */
	ya = y[j];
    }
} /* seife_rekfl */


/* ----- END OF cseife_rekfl.c ----- */
