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

/*------------------------------------------------------------------------
 * various interpolation routines;
 * the sinc and spline functions are not thread-safe
 *
 * sinc interpolation adapted from Seismic Unix, CWP
 * spline interpolation adapted from Numerical Recipes in C
 *------------------------------------------------------------------------*/

#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>

#include "interpol.h"

#define NTABLE 513
#define LTABLE 8

static const float fntablem1 = 512.0f;
static bool b_sinc_initialized = false;

static float table[NTABLE][LTABLE];

static inline float sincf(float x)
{
  if (0.0 == x) {
    return 1.0;
  } else {
    float pix = M_PI * x;
    return sinf(pix)/pix;
  }
}

static void solve(float *r, float *g, float *f, float *a)
{
  float v, e, c, w, bot;

  if (r[0] == 0.0) return;

  a[0] = 1.0f;
  v = r[0];
  f[0] = g[0]/r[0];

  for (int j=1; j<LTABLE; ++j) {
    // solve Ra=v, see Claerbout, FGDP, p.57
    a[j] = 0.0f;
    f[j] = 0.0f;
    e = 0.0f;
    for (int i=0; i<j; ++i) e += a[i]*r[j-i];
    c = e/v;
    v -= c*e;
    for (int i=0; i<=j/2; ++i) {
      bot = a[j-i]-c*a[i];
      a[i] -= c*a[j-i];
      a[j-i] = bot;
    }
    w = 0.0f;
    for (int i=0; i<j; ++i) w += f[i]*r[j-i];
    c = (w-g[j])/v;
    for (int i=0; i<=j; ++i) f[i] -= c*a[j-i];
  }

  return;
}

static void mksinc(float d, float *sinc)
{
  static float ss[20], aa[20], cc[20], ww[20];

  float *s = &ss[0];
  float *a = &aa[0];
  float *c = &cc[0];
  float *work = &ww[0];
  float fmax;

  // compute auto-correlation and cross-correlation arrays
  fmax = 0.066+0.265*logf((float)LTABLE);
  fmax = (fmax<1.0f) ? fmax : 1.0f;
  for (int j=0; j<LTABLE; ++j) {
    a[j] = sincf(fmax*j);
    c[j] = sincf(fmax*(LTABLE/2-j-1+d));
  }
  // solve symmetric Toeplitz system for the sinc approximation
  solve(a,c,s,work);
  for (int j=0; j<LTABLE; ++j) sinc[j] = s[j];

  return;
}

void sinc_init()
{
  float frac, *p;

  for (int jtable=1; jtable<NTABLE-1; ++jtable) {
    frac = (float)jtable/fntablem1;
    p = &table[jtable][0];
    mksinc(frac,p);
  }
  for (int jtable=0; jtable<LTABLE; ++jtable) {
    table[0][jtable] = 0.0f;
    table[NTABLE-1][jtable] = 0.0f;
  }
  table[0][LTABLE/2-1] = 1.0f;
  table[NTABLE-1][LTABLE/2] = 1.0f;

  b_sinc_initialized = true;

  return;
}

void sinc_ipol(int ntin,        // input number of samples at which ampin is given
	       float dtin,      // input sampling interval dt
	       float ftin,      // first input time sample (usually ftin = 0.0)
	       float *ampin,    // input amplitude values; note: ampin[0] = amp[ftin]
	       int ntout,       // output number of samples
	       float *tout,     // time values at which ampout is calculated
	       float *ampout)   // output amplitude values; note: ampout[0] = amp[tout[0]]
{
  if (!b_sinc_initialized) {
    sinc_init();
  }

  int ioutb, ntinm8, itoutn, kampin, ktable;
  float toutb, toutf, touts, toutn, frac, ampini;
  float sum, *ampin0, *table00, *pampin, *ptable;

  ioutb = -11;
  toutf = ftin;
  touts = 1.0/dtin;
  toutb = 8.0-toutf*touts;
  ntinm8 = ntin-8;
  ampin0 = &ampin[0];
  table00 = &table[0][0];

  // loop over output samples
  for (int itout=0; itout<ntout; ++itout) {
    // determine pointers into table and yin
    toutn = toutb+tout[itout]*touts;
    itoutn = (int)toutn;
    kampin = ioutb+itoutn;
    pampin = ampin0+kampin;
    frac = toutn-(float)itoutn;
    ktable = (frac>=0.0) ? (int)(frac*fntablem1+0.5) : (int)((frac+1.0)*fntablem1-0.5);
    ptable = table00+ktable*8;

    // if totally within input array, use fast method
    if (kampin>=0 && kampin<=ntinm8) {
      ampout[itout] =
        pampin[0]*ptable[0]+
        pampin[1]*ptable[1]+
        pampin[2]*ptable[2]+
        pampin[3]*ptable[3]+
        pampin[4]*ptable[4]+
        pampin[5]*ptable[5]+
        pampin[6]*ptable[6]+
        pampin[7]*ptable[7];
      // else handle end effects with care
    } else {
      // sum over 8 tabulated coefficients
      sum = 0.0f;
      for (int itable=0; itable<8; ++itable, ++kampin) {
	// mirror trace at the beginning while negating amplitude
	if (kampin<0) ampini = -ampin[-kampin];
	// don't mirror at the end as our signature should be zero (close to zero) anyway
	else if (kampin>=ntin) ampini = 0.0f; 
        else ampini = ampin[kampin];
        sum += ampini*(*ptable++);
      }
      ampout[itout] = sum;
    }
  }

  return;
}

void lin_ipol(int ntin,        // input number of samples at which ampin is given
	      float dtin,      // input sampling interval dt
	      float ftin,      // first input time sample (usually ftin = 0.0)
	      float *ampin,    // input amplitude values; note: ampin[0] = amp[ftin]
	      int ntout,       // output number of samples
	      float *tout,     // time values at which ampout is calculated
	      float *ampout)   // output amplitude values; note: ampout[0] = amp[tout[0]]
{
  float idxinf, w;
  int idxin;

  for (int itout=0; itout<ntout; ++itout) {
    idxinf = (tout[itout]-ftin)/dtin;
    idxin = floorf(idxinf);
    if (idxin < 0) {
      ampout[itout] = ampin[0];
    } else if (idxin >= ntin-1) {
      ampout[itout] = ampin[ntin-1];
    } else {
      w = idxinf - (float)idxin;
      ampout[itout] = w*ampin[idxin+1] + (1.0f-w)*ampin[idxin];
    } 
  }

  return;
}

void spline_init(float *x,   // x values of points known
		 float *y,   // y values of points known
		 int n,      // number of x-y pairs
		 float yp1,  // 1st derivative at front (>= 10^30 to do nat'l spline)
		 float ypn,  // 1st derivative at back (>= 10^30 to do nat'l spline)
		 float *y2)  // output, 2nd derivative to use in splint
{
  int k;
  float	p, qn, sig, un;

  float *u = (float*)malloc((n-1)*sizeof(float));

  if(yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else{
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for(int i = 1; i < n-1; i++){
    sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
  }
  if(ypn > 0.99e30)
    qn = un = 0.0;
  else{
    qn = 0.5;
    un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
  }
  y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);
  for(k = n-2; k >= 0; k--){
    y2[k] = y2[k]*y2[k+1] + u[k];
  }

  free(u);

  return;
}

float spl_ipol(float *xa,    // x values of points known
	       float *ya,    // y values of points known
	       float *y2a,   // 2nd derivative from spline
	       int n,        // number of x-y pairs
	       float x)      // x value where to interpolate
{
  int klo, khi, k;
  float	h, b, a;
  static int pklo = 0, pkhi = 1;

  if(xa[pklo] <= x && xa[pkhi] > x) {
    klo = pklo;
    khi = pkhi;
  } else {
    klo = 0;
    khi = n - 1;
    while(khi - klo > 1) {
      k = (khi + klo) >> 1;
      if(xa[k] > x) khi = k;
      else          klo = k;
    }
    pklo = klo;
    pkhi = khi;
  }

  h = xa[khi] - xa[klo];
  assert(h != 0); // otherwise user supplied rubbish xa data with duplicate points
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;

  return a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}
