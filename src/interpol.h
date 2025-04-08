/*------------------------------------------------------------------------
 * Copyright (C) 2024 For the list of authors, see file AUTHORS.
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
 * along with IFOS. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 *------------------------------------------------------------------------*/

#ifndef _INTERPOL_H_
#define _INTERPOL_H_

// initialize sinc interpolation; if not called explicitly up front, it
// will be called at the beginning of the first sinc_ipol function call;
// the sinc functions are not thread-safe
void sinc_init();

// sinc interpolation
void sinc_ipol(int ntin,        // input number of samples at which ampin is given
	       float dtin,      // input sampling interval dt
	       float ftin,      // first input time sample (usually ftin = 0.0)
	       float *ampin,    // input amplitude values; note: ampin[0] = amp[ftin]
	       int ntout,       // output number of samples
	       float *tout,     // time values at which ampout is calculated
	       float *ampout);  // output amplitude values; note: ampout[0] = amp[tout[0]]

// linear interpolation
void lin_ipol(int ntin,        // input number of samples at which ampin is given
	      float dtin,      // input sampling interval dt
	      float ftin,      // first input time sample (usually ftin = 0.0)
	      float *ampin,    // input amplitude values; note: ampin[0] = amp[ftin]
	      int ntout,       // output number of samples
	      float *tout,     // time values at which ampout is calculated
	      float *ampout);  // output amplitude values; note: ampout[0] = amp[tout[0]]

// initialize spline interpolation; needs to be called for each trace to
// be interpolated; the spline functions are not thread-safe
void spline_init(float *x,    // x values of points known
		 float *y,    // y values of points known
		 int n,       // number of x-y pairs
		 float yp1,   // 1st derivative at front (>= 10^30 to do nat'l spline)
		 float ypn,   // 1st derivative at back (>= 10^30 to do nat'l spline)
		 float *y2);  // output, 2nd derivative to use in splint

// spline interpolation; return value: interpolated value y(x)
float spl_ipol(float *xa,      // x values of points known
	       float *ya,      // y values of points known
	       float *y2a,     // 2nd derivative from spline
	       int n,          // number of x-y pairs
	       float x);       // x value where to interpolate

#endif
