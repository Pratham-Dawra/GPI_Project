
/*------------------------------------------------------------------------
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
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating particle velocities at interior gridpoints (excluding boundarys) [GX2+1...GX3][GY2+1...GY3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_v_interior(int nt, float **srcpos_loc, float **signals, int nsrc,
                       MemModel *mpm, MemWavefield *mpw, GlobVar *gv)
{
    float amp;
    float sxx_x, sxy_x, sxy_y, syy_y;
    double time1 = 0.0, time2 = 0.0;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating particle velocities...\n");
    }

    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */

    for (int l = 1; l <= nsrc; l++) {
        int i = (int)srcpos_loc[1][l];
        int j = (int)srcpos_loc[2][l];
        float azi_rad = srcpos_loc[7][l] * PI / 180;

        //amp=signals[l][nt]; // unscaled force amplitude
        amp = (gv->DT * signals[l][nt]) / (gv->DH * gv->DH);    // scaled force amplitude with F= 1N

        gv->SOURCE_TYPE = (int)srcpos_loc[8][l];

        switch (gv->SOURCE_TYPE) {
          case 2:              /* single force in x */
              mpw->pvx[j][i] += mpm->prip[j][i] * amp;
              /* previous implementation of body forces as seismic sources.
               * Implementation according to Coutant et al., BSSA, Vol. 85, No 5, 1507-1512.
               * The stress tensor components sxx and syy are incremented prior to
               * particle velocity update. Thereby the body force (both directions)
               * are located at full grid point (i,j) (same position as pressure source).
               * as a consequence, source signals are added [and weighted] at multiple grid points.
               * This implementation works but not quite physical when considering e.g. a force of 1 N
               * and aiming to gain the particle velocity strictly according to that force.
               * Therefore it has been commented */

              /*for (m=1; m<=fdoh; m++) {
               * vx[j][i+m-1]  +=  hc[m]*rip[j][i]*amp;
               * vx[j][i-m]    +=  hc[m]*rip[j][i-1]*amp;
               * } */
              break;
          case 3:              /* single force in y */
              mpw->pvy[j][i] += mpm->prjp[j][i] * amp;
              /*for (m=1; m<=fdoh; m++) {
               * vy[j+m-1][i]  +=  hc[m]*rjp[j][i]*amp;
               * vy[j-m][i]    +=  hc[m]*rjp[j][i-1]*amp;
               * } */
              break;
          case 4:              /* custom force */
              mpw->pvx[j][i] += sin(azi_rad) * (mpm->prip[j][i] * amp);
              mpw->pvy[j][i] += cos(azi_rad) * (mpm->prjp[j][i] * amp);
              /*for (m=1; m<=fdoh; m++) {
               * vx[j][i+m-1]  +=  sin(azi_rad)*(hc[m]*rip[j][i]*amp);
               * vx[j][i-m]    +=  sin(azi_rad)*(hc[m]*rip[j][i-1]*amp);
               * vy[j+m-1][i]  +=  cos(azi_rad)*(hc[m]*rjp[j][i]*amp);
               * vy[j-m][i]    +=  cos(azi_rad)*(hc[m]*rjp[j][i-1]*amp);
               * } */
              break;
        }
    }

  
    
    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {  /* elastic cases */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                gv->FDOP_V(i, j, &sxx_x, &sxy_x, &sxy_y, &syy_y, mpw);
                wavefield_update_v(i, j, sxx_x, sxy_x, sxy_y, syy_y, mpm, mpw, gv);
            }
        }
        
    }else /*acoustic cases */
    {
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                gv->FDOP_AC_V(i, j, &sxx_x, &syy_y, mpw);
                wavefield_update_v_ac(i, j, sxx_x, syy_y, mpm, mpw, gv);
            }
        }
    }

    
    
    
    
    
    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating particle velocities (real time: %4.3fs).\n", time2 - time1);
    }
}
