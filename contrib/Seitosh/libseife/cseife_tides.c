/*! \file cseife_tides.c
 * \brief remove tides (implementation)
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
 * \date 22/06/2020
 * 
 * remove tides (implementation)
 * 
 * REVISIONS and CHANGES 
 *  - 28/06/2005   V1.0   Thomas Forbriger
 *  - 22/06/2020   V1.1   handle case of nstep=1 explicitly
 *                 V1.2   make array size flexible and add M3 frequency
 * 
 * ============================================================================
 */
#define TF_CSEIFE_TIDES_C_VERSION \
  "TF_CSEIFE_TIDES_C   V1.1"

#include "cseife.h"
#include <math.h>

/* subs/seife_tides.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

   the code was derived through f2c, but modified thereafter
*/

/* detide with synthetic tides interpolated over ni samples */
void seife_tid(double* x, int n, double dt, int nstep)
{
    /* Initialized data */

/* tidal constituents used:
 *
 * f /cpd         f / deg/h   constituent
 * ------         ---------   -----------
 * 0.89293        13.39395	Q1
 * 0.92954        13.9431	O1
 * 1.00274        15.0411	K1
 * 1.89567        28.43505	N2
 * 1.93227        28.98405	M2
 * 2              30          S2
 * 2.898450424    43.47675636 M3
*/

#define C_MFREQ 7
#define C_MDIM (C_MFREQ*2+1)
    static double omega[C_MFREQ] = { 2.898450424,
      1.93227,.92954,2.,1.00274,1.89567,.89293 };
    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Builtin functions */
    // double atan(double), cos(double), sin(double);

    /* Local variables */
    static double tdif, omeg[C_MFREQ];
    static int ndim;
    static double step, tint, omeg0, xgez1, xgez2, xgez3, 
                  a[C_MDIM*C_MDIM], c__[C_MDIM], d__[C_MDIM], 
                  e[C_MDIM], f[C_MDIM];
    static int i__, j, k;
    static double t;
    static int nfreq, k2, k3;
    static double tstep2;
    static int jj;
    static double rs[C_MDIM], sx, dth;
    static int nco;
    static double cor[C_MFREQ], dur;
//    extern /* Subroutine */ int seife_gauss__(doublereal *, integer *, 
//	    integer *, doublereal *, doublereal *);

/*  remove tides. number of frequencies is automatically chosen according */
/*  to the total length of the record. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ndim = C_MDIM;
    dth = dt / two;

    if (nstep == 0) {
	nstep = (float)300. / dt;
    }
    nstep = nstep > 1 ? nstep : 1;
    step = (double) nstep;
    tstep2 = (step - one) * dth;
    tint = tstep2 + dth;
    if (nstep == 1) {
	step = 1.;
	tstep2 = 0.;
	tint = dth;
    }
/*  determine the number of frequencies required for a good fit */
    dur = n * dt / 3600.;
    nfreq = C_MFREQ;
    if (dur < 35.) {
	nfreq = 6;
    }
    if (dur < 18.) {
	nfreq = 5;
    }
    if (dur < 14.) {
	nfreq = 4;
    }
    if (dur < 5.) {
	nfreq = 3;
    }
    nco = (nfreq << 1) + 1;
    omeg0 = atan(one) * 8. / 86400.;
    i__1 = nfreq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L101: */
	omeg[i__ - 1] = omeg0 * omega[i__ - 1];
    }
/*  determine partial amplitudes by least-squares fit */
    i__1 = nco;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rs[i__ - 1] = 0.;
	i__2 = nco;
	for (k = 1; k <= i__2; ++k) {
/* L1: */
	    a[i__ + k * C_MDIM - (C_MDIM+1)] = 0.;
	}
    }
    if (nstep == 1) {
/* ====================================================================== */
/* do not average */
/* ============== */
/*  correction for averaging over nstep samples */
	i__2 = nfreq;
	for (j = 1; j <= i__2; ++j) {
	    cor[j - 1] = 1.;
	}
/*  set up system of linear equations */
	c__[0] = one;
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    sx = zero;
	    i__1 = j;
	    for (jj = j; jj <= i__1; ++jj) {
		sx += x[jj];
	    }
	    t = (j - 1) * dt;
	    i__1 = nfreq;
	    for (k = 1; k <= i__1; ++k) {
		c__[(k << 1) - 1] = cos(omeg[k - 1] * t);
		c__[k * 2] = sin(omeg[k - 1] * t);
	    }
	    i__1 = nco;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		rs[i__ - 1] += sx * c__[i__ - 1];
		i__3 = nco;
		for (k = 1; k <= i__3; ++k) {
		    a[i__ + k * C_MDIM - (C_MDIM+1)] += c__[i__ - 1] * c__[k - 1];
		}
	    }
	}
/*  solve for partial amplitudes */
	seife_gauss(a, nco, ndim, rs, f);
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    t = (j - 1) * dt;
/*  remove average and tides */
	    xgez1 = f[0];
	    i__1 = nfreq;
	    for (k = 1; k <= i__1; ++k) {
		k2 = k << 1;
		k3 = k2 + 1;
		c__[k2 - 1] = cos(omeg[k - 1] * t);
		c__[k3 - 1] = sin(omeg[k - 1] * t);
		xgez1 = xgez1 + f[k2 - 1] * c__[k2 - 1] + f[k3 - 1] * c__[k3 
			- 1];
	    }
	    x[j] -= xgez1;
	}
    } else {
/* ====================================================================== */
/* average over nstep samples */
/* ========================== */
/*  correction for averaging over nstep samples */
	i__2 = nfreq;
	for (j = 1; j <= i__2; ++j) {
/* L102: */
	    cor[j - 1] = step * sin(omeg[j - 1] * dth) / sin(step * omeg[j - 
		    1] * dth);
	}
/*  set up system of linear equations */
	c__[0] = one;
	i__2 = n - nstep + 1;
	i__1 = nstep;
	for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	    sx = zero;
	    i__3 = j + nstep - 1;
	    for (jj = j; jj <= i__3; ++jj) {
/* L12: */
		sx += x[jj];
	    }
	    sx /= step;
	    t = (j - 1) * dt + tstep2;
	    i__3 = nfreq;
	    for (k = 1; k <= i__3; ++k) {
		c__[(k << 1) - 1] = cos(omeg[k - 1] * t);
/* L103: */
		c__[k * 2] = sin(omeg[k - 1] * t);
	    }
	    i__3 = nco;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		rs[i__ - 1] += sx * c__[i__ - 1];
		i__4 = nco;
		for (k = 1; k <= i__4; ++k) {
/* L2: */
		    a[i__ + k * C_MDIM - (C_MDIM+1)] += c__[i__ - 1] * c__[k - 1];
		}
	    }
	}
/*  solve for partial amplitudes */
	seife_gauss(a, nco, ndim, rs, f);
	i__4 = n - nstep + 1;
	i__3 = nstep;
	for (j = 1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
	    t = (j - 1) * dt + tstep2;
/*  remove average and tides */
	    i__1 = nfreq;
	    for (k = 1; k <= i__1; ++k) {
		k2 = k << 1;
		k3 = k2 + 1;
		c__[k2 - 1] = cos(omeg[k - 1] * t);
		c__[k3 - 1] = sin(omeg[k - 1] * t);
		d__[k2 - 1] = -omeg[k - 1] * c__[k3 - 1];
		d__[k3 - 1] = omeg[k - 1] * c__[k2 - 1];
		e[k2 - 1] = -omeg[k - 1] * d__[k3 - 1];
/* L104: */
		e[k3 - 1] = omeg[k - 1] * d__[k2 - 1];
	    }
	    xgez1 = f[0];
	    xgez2 = zero;
	    xgez3 = zero;
	    i__1 = nfreq;
	    for (k = 1; k <= i__1; ++k) {
		k2 = k << 1;
		k3 = k2 + 1;
		xgez1 += cor[k - 1] * (f[k2 - 1] * c__[k2 - 1] + f[k3 - 1] * 
			c__[k3 - 1]);
		xgez2 += cor[k - 1] * (f[k2 - 1] * d__[k2 - 1] + f[k3 - 1] * 
			d__[k3 - 1]);
/* L105: */
		xgez3 += cor[k - 1] * (f[k2 - 1] * e[k2 - 1] + f[k3 - 1] * e[
			k3 - 1]) / two;
	    }
	    if (j > n - (nstep << 1) + 1) {
		nstep = n + 1 - j;
	    }
	    i__1 = j + nstep - 1;
	    for (jj = j; jj <= i__1; ++jj) {
		tdif = (jj - j) * dt - tstep2;
/* L3: */
		x[jj] = x[jj] - xgez1 - xgez2 * tdif - xgez3 * tdif * tdif;
	    }
	}
    }
} /* seife_tides__ */


/* ----- END OF cseife_tides.c ----- */
