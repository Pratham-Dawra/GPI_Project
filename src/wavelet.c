
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

/*------------------------------------------------------------------------
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_srcsig.h"
#include <complex.h>

float **wavelet(float **srcpos_loc, int nsrc, GlobVar *gv)
{
    int nts = 0;
    float *psource = NULL, tshift, amp = 0.0, a, fc, tau, t, ts, f2;
    float n, alpha, phi0deg, phi0, T, fmin, fmax, width, taper, sigmax;
    float **signals = NULL;

    if (gv->SOURCE_SHAPE == 3)
        psource = read_srcsig(&nts, gv);

    /* If there is no source in the current domain, return early as otherwise
     * the call to matrix with nsrc=0 will cause memory corruption. It would
     * determine nrow=0 and allocate one slot (C-index 0), later however try
     * to access C-index 1 to start allocating NT columns. */
    if (nsrc < 1)
        return signals;

    signals = matrix(1, nsrc, 1, gv->NT);

    for (int k = 1; k <= nsrc; k++) {
        tshift = srcpos_loc[4][k];  // time shift
        fc = srcpos_loc[5][k];  // centre frequency [Hz]; in case 5 (Berlage wavelet) lowest frequency
        a = srcpos_loc[6][k];   // maximum amplitude
        ts = 1.0 / fc;
        sigmax = 0.0;

        for (int nt = 1; nt <= gv->NT; nt++) {
            t = (float)nt *gv->DT;

            switch (gv->SOURCE_SHAPE) {
              case 1:
                  /* standard Ricker signal */
                  tau = PI * (t - 1.5 * ts - tshift) / (ts);
                  amp = (((1.0 - 2.0 * tau * tau) * exp(-tau * tau)));
                  break;
              case 2:
                  if ((t < tshift) || (t > (tshift + ts)))
                      amp = 0.0;
                  else
                      amp = (sin(2.0 * PI * (t - tshift) * fc) - 0.5 * sin(4.0 * PI * (t - tshift) * fc));
                  break;
              case 3:
                  if (nt <= nts)
                      amp = psource[nt - 1];    // psource has zero-based index
                  else
                      amp = 0.0;
                  break;        /* source wavelet from file SOURCE_FILE */
              case 4:
                  if ((t < tshift) || (t > (tshift + ts)))
                      amp = 0.0;
                  else
                      amp = pow(sin(PI * (t - tshift) / ts), 3.0);
                  break;        /* sinus raised to the power of three */
              case 5:
                  /* Berlage wavelet (minimum-phase) (Aldridge, 1990) */
                  n = srcpos_loc[9][k];         // time exponent; >0 (Berlage only)
                  alpha = srcpos_loc[10][k];    // exponential decay factor (Berlage only)
                  phi0deg = srcpos_loc[11][k];  // initial phase angle [°] (Berlage only)
                  phi0 = phi0deg * PI / 180;
                  if (n <= 0)
                      log_fatal("No valid time exponent for Berlage wavelet (N>0) specified!\n");
                  if ((t < tshift))
                      amp = 0.0;
                  else
                      amp =
                          a * pow((t - tshift),
                                  n) * exp(-alpha * (t - tshift)) * cos(2 * PI * fc * (t - tshift) + phi0);
                  break;
              case 6:
                  /* Klauder wavelet */
                  fmin = srcpos_loc[9][k];   // lowest frequency in sweep [Hz] (Klauder only)
                  fmax = srcpos_loc[10][k];  // highest frequency in sweep [Hz] (Klauder only)
                  T = srcpos_loc[11][k];     // sweep duration [s] (Klauder only)
                  width = srcpos_loc[12][k]; // width of the Klauder wavelet in number of centre periods

                  f2 = (fmax - fmin) / T;
                  tau = (width / fc);   // Klauder wavelet is shifted by "width" centre periods
                  if ((t < tshift) || (t > (tau * 2) + tshift)) {
                      amp = 0.0;
                  } else if (t < (tau * 0.2) + tshift) {
                      amp =
                          creal(sin(PI * f2 * (t - tau - tshift) * (T - (t - tau - tshift))) /
                                (PI * f2 * (t - tau - tshift)) * exp(2 * PI * I * fc * (t - tau - tshift)));
                      taper = (1 - cos(((t - tshift) / tau / 0.2) * PI)) / 2;   /* linear taper */
                      amp = amp * taper;
                  } else if (t > (tau * 1.8) + tshift) {
                      amp =
                          creal(sin(PI * f2 * (t - tau - tshift) * (T - (t - tau - tshift))) /
                                (PI * f2 * (t - tau - tshift)) * exp(2 * PI * I * fc * (t - tau - tshift)));
                      taper = (1 - cos(((2 * tau + tshift - t) / 0.2 / tau) * PI)) / 2; /* linear taper */
                      amp = amp * taper;
                  } else {
                      amp =
                          creal(sin(PI * f2 * (t - tau - tshift) * (T - (t - tau - tshift))) /
                                (PI * f2 * (t - tau - tshift)) * exp(2 * PI * I * fc * (t - tau - tshift)));
                  }
                  break;
              default:
                  log_fatal("No valid source wavelet (gv->SOURCE_SHAPE) specified!\n");
                  break;
            }
            if (fabs(amp) > sigmax)
                sigmax = fabs(amp);
            signals[k][nt] = amp * a;
        }
        /* Normalization to desired amplitude */
        if (0.0 == sigmax)
            log_fatal("Source signal contains only zeros!\n");
        if (gv->SOURCE_SHAPE == 5 || gv->SOURCE_SHAPE == 6) {
            for (int nt = 1; nt <= gv->NT; nt++) {
                signals[k][nt] = signals[k][nt] / sigmax;
            }
        }
    }

    log_info("%d source position(s) in subdomain assigned with source signal.\n", nsrc);

    if (gv->SOURCE_SHAPE == 3)
        free(psource);

    return signals;
}
