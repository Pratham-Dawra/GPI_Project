
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

/*-------------------------------------------------------------
 *  Check FD-Grid for stability and grid dispersion.
 *  If the stability criterion is not fulfilled the program will
 *  terminate.
 *
 *  ----------------------------------------------------------*/

#include <limits.h>

#include "fd.h"
#include "logging.h"

/*************************************************************
 * TODO: Calculation of global grid size wrong! To be fixed. *
 *************************************************************/

void checkfd(float *hc, float **srcpos, int nsrc, int **recpos, GlobVar * gv)
{
    float fmax, gamma;
    float dtstab, dhstab, temporal;
    float snapoutx = 0.0, snapouty = 0.0;
    float srec_minx = gv->DH * gv->NX * gv->NPROCX + 1, srec_miny = gv->DH * gv->NY * gv->NPROCY + 1;
    float srec_maxx = -1.0, srec_maxy = -1.0;
    float CFL;

    if (0 == gv->MPID) {
        log_info("------------------------- Min/max velocities ----------------\n");
        log_info("Note: If any velocity is set <1m/s, it will be ignored while\n");
        log_info("determining stable DH and DT parameters.\n");
    }

    float cmax = gv->VSMAX > gv->VPMAX ? gv->VSMAX : gv->VPMAX;
    float cmin = gv->VSMIN < gv->VPMIN ? gv->VSMIN : gv->VPMIN;

    if (gv->FDORDER_TIME == 4) {
        temporal = 3.0 / 2.0;
    } else {
        temporal = 1.0;
    }

    fmax = 2.0 / gv->TS;
    dhstab = (cmin / (hc[0] * fmax));
    gamma = fabs(hc[1]) + fabs(hc[2]) + fabs(hc[3]) + fabs(hc[4]) + fabs(hc[5]) + fabs(hc[6]);
    dtstab = gv->DH / (sqrt(2) * gamma * cmax * temporal);
    CFL = cmax * gv->DT / gv->DH;

    if (gv->MPID == 0) {
        log_info("The following velocities take anisotropy and/or attenuation into account.\n");
        log_info("Min and max P-wave phase velocity: Vp_min=%.2fm/s, Vp_max=%.2fm/s\n", gv->VPMIN, gv->VPMAX);
        log_info("Min and max S-wave phase velocity: Vs_min=%.2fm/s, Vs_max=%.2fm/s\n", gv->VSMIN, gv->VSMAX);
        log_info("Overall global min and max velocity: V_min=%.2fm/s, V_max=%.2fm/s\n", cmin, cmax);
        log_info("------------------------- Grid dispersion -------------------\n");
        log_info("To limit grid dispersion, the number of grid points per min.\n");
        log_info("wavelength (of S-waves) should be at least 6. Here, the min.\n");
        log_info("wavelength is assumed to be the minimum model phase velocity\n");
        log_info("(of S-waves) at max. frequency divided by the max. frequency\n");
        log_info("of the source (here approximately %.3fHz). The minimum\n", 2.0 / gv->TS);
        log_info("wavelength in the following simulation will be %.3fm. Thus,\n", cmin / fmax);
        log_info("the recommended value is DH=%.3fm. Your value: DH=%.3fm\n", dhstab, gv->DH);
        if (gv->DH > dhstab) {
            log_warn("Grid dispersion will influence wave propagation, choose a smaller grid spacing (DH).\n");
        }

        log_info("------------------------- Stability check -------------------\n");
        log_info("The simulation is stable if p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n");
        log_info("where cmax is the maximum phase velocity at infinite frequency\n");
        log_info("and gamma=sum(|FD coeff.|). In this simulation, cmax=%.2fm/s.\n", cmax);
        log_info("DT is the time step and DH is the grid size. The Courant-\n");
        log_info("Friedrichs-Lewy number is %.4f. The stability limit for the\n", CFL);
        log_info("time step is DT=%es. Your value: DT=%es\n", dtstab, gv->DT);
        if (gv->DT > dtstab) {
            log_error("The simulation will get unstable, choose a smaller DT.\n");
            log_fatal("Instability of the simulation, update your parameters.\n");
        } else {
            log_info("The simulation will be stable.\n");
        }

        log_info("------------------------- Additional checks -----------------\n");
        if (gv->SNAP) {
            if (gv->TSNAP2 > gv->TIME) {
                log_warn
                    ("TSNAP2=%e (last snapshot) > time of wave propagation %e. TSNAP2 changed to be equal to TIME.\n",
                     gv->TSNAP2, gv->TIME);
                gv->TSNAP2 = gv->TIME;
            }
            snapoutx = gv->NX / (float)gv->IDX;
            snapouty = gv->NY / (float)gv->IDY;
            log_info("Output of snapshot grid points (x) per node: %8.2f\n", snapoutx);
            log_info("Output of snapshot grid points (y) per node: %8.2f\n", snapouty);
            if (snapoutx - (int)snapoutx > 0)
                log_fatal("Ratio NX-NPROCX-IDX must be integer.\n");
            if (snapouty - (int)snapouty > 0)
                log_fatal("Ratio NY-NPROCY-IDY must be integer.\n");
        } else {
            log_info("Skipping checks of snapshot parameters.\n");
        }

        if (gv->SEISMO) {
            log_info("Number of modeling time steps: %d\n", gv->NT);
            log_info("Seismogram sampling interval in time steps: %d\n", gv->NDT);
            log_info("Number of seismogram output samples: %d\n", gv->NT / gv->NDT);

            /* SU and SEG-Y allow 32767 samples, furthermore the exist programs allow for 65535 
             * samples and pseudo-SEG-Y formats allow foralmost arbitrarily long traces.
             * For binary and textual output the limit is arbitrary. USHRT_MAX is the limut of 
             * an unsigned short specified in limits.h */

            if ((gv->SEIS_FORMAT == 1) && (gv->NT / gv->NDT) > (USHRT_MAX)) {
                log_error("Maximum number of samples per trace in SU format: %d. Your value: %d\n", USHRT_MAX,
                          gv->NT / gv->NDT);
                log_fatal("Too many output samples per receiver for SU format.\n");
            }

            srec_minx = gv->DH * gv->NX * gv->NPROCX + 1, srec_miny = gv->DH * gv->NY * gv->NPROCY + 1;
            srec_maxx = -1.0, srec_maxy = -1.0;
            log_info("Checking for receiver position(s) as specified in json file.\n");
            log_info("Global grid size in m: %5.2f (x) : %5.2f (y)\n", gv->NX * gv->DH * gv->NPROCX,
                     gv->NY * gv->DH * gv->NPROCY);
            if (gv->FREE_SURF == 0)
                log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n",
                         (float)gv->FW * gv->DH, gv->NX * gv->DH * gv->NPROCX - (float)gv->FW * gv->DH,
                         (float)gv->FW * gv->DH, gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH);
            if (gv->FREE_SURF == 1)
                log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n",
                         (float)gv->FW * gv->DH, gv->NX * gv->DH * gv->NPROCX - (float)gv->FW * gv->DH, gv->DH,
                         gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH);

            /* find maximum and minimum source positions coordinate ---- from input file */
            if (gv->READREC == 0) {
                if (gv->XREC1 > gv->XREC2) {
                    srec_maxx = gv->XREC1;
                    srec_minx = gv->XREC2;
                } else {
                    srec_maxx = gv->XREC2;
                    srec_minx = gv->XREC1;
                }
                if (gv->YREC1 > gv->YREC2) {
                    srec_maxy = gv->YREC1;
                    srec_miny = gv->YREC2;
                } else {
                    srec_maxy = gv->YREC2;
                    srec_miny = gv->YREC1;
                }
                log_info("Number of receiver positions: %d\n", gv->NTRG);
            }

            if (gv->READREC == 1) {
                /* find maximum and minimum source positions coordinate ---- from receiver file */
                for (int k = 1; k <= gv->NTRG; k++) {
                    /* find maximum source positions coordinate */
                    if ((recpos[1][k] * gv->DH) > srec_maxx)
                        srec_maxx = recpos[1][k] * gv->DH;
                    if ((recpos[2][k] * gv->DH) > srec_maxy)
                        srec_maxy = recpos[2][k] * gv->DH;
                    /* find minimum source positions coordinate */
                    if ((recpos[1][k] * gv->DH) < srec_minx)
                        srec_minx = recpos[1][k] * gv->DH;
                    if ((recpos[2][k] * gv->DH) < srec_miny)
                        srec_miny = recpos[2][k] * gv->DH;
                }
                log_info("Number of receiver positions: %d\n", gv->NTRG);
            }

            log_info("Minimum receiver position coordinates: %5.2f (x) : %5.2f (y)\n", srec_minx, srec_miny);
            log_info("Maximum receiver position coordinates: %5.2f (x) : %5.2f (y)\n", srec_maxx, srec_maxy);

            /* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
            if (((srec_maxx < 0.0) || (srec_maxy < 0.0)) || ((srec_minx < 0.0) || (srec_miny < 0.0))) {
                log_fatal("Coordinate of at least one receiver location is outside the global grid.\n");
            }
            if ((srec_maxx > gv->NX * gv->DH * gv->NPROCX) || (srec_maxy > gv->NY * gv->DH * gv->NPROCY)) {
                log_fatal("Coordinate of at least one receiver location is outside the global grid.\n");
            }

            /* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
            if ((srec_maxx < ((float)gv->FW * gv->DH)) || (srec_minx < ((float)gv->FW * gv->DH))) {
                /* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary") */
                log_warn
                    ("Coordinate of at least one receiver location is inside the Absorbing Boundary (left boundary).\n");
            }

            if (srec_maxx > (gv->NX * gv->DH * gv->NPROCX - (float)gv->FW * gv->DH)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one receiver location is inside the Absorbing Boundary (right boundary).\n");
            }

            if (srec_maxy > (gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one receiver location is inside the Absorbing Boundary (lower boundary).\n");
            }

            if ((srec_miny < ((float)gv->FW * gv->DH)) && !(gv->FREE_SURF)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one receiver location is inside the Absorbing Boundary (top boundary).\n");
            }
        } else {
            log_info("Skipping checks of seismogram parameters.\n");
        }

        if (gv->SRCREC == 1) {
            srec_minx = gv->DH * gv->NX * gv->NPROCX + 1, srec_miny = gv->DH * gv->NY * gv->NPROCY + 1;
            srec_maxx = -1.0, srec_maxy = -1.0;
            log_info("Checking for source position(s) as specified in source file.\n");
            log_info("Global grid size in m: %5.2f (x) : %5.2f (y)\n", gv->NX * gv->DH * gv->NPROCX,
                     gv->NY * gv->DH * gv->NPROCY);
            if (gv->FREE_SURF == 0)
                log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n",
                         (float)gv->FW * gv->DH, gv->NX * gv->DH * gv->NPROCX - (float)gv->FW * gv->DH,
                         (float)gv->FW * gv->DH, gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH);
            if (gv->FREE_SURF == 1)
                log_info("Global grid size in m(-width of abs.boundary): %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m)\n",
                         (float)gv->FW * gv->DH, gv->NX * (float)gv->DH * gv->NPROCX - (float)gv->FW * gv->DH, gv->DH,
                         gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH);

            for (int k = 1; k <= nsrc; k++) {
                /* find maximum source positions coordinate */
                if (srcpos[1][k] > srec_maxx)
                    srec_maxx = srcpos[1][k];
                if (srcpos[2][k] > srec_maxy)
                    srec_maxy = srcpos[2][k];
                /* find minimum source positions coordinate */
                if (srcpos[1][k] < srec_minx)
                    srec_minx = srcpos[1][k];
                if (srcpos[2][k] < srec_miny)
                    srec_miny = srcpos[2][k];
            }

            log_info("Number of source positions: %d\n", nsrc);
            log_info("Minimum source position coordinates: %5.2f (x) : %5.2f (y)\n", srec_minx, srec_miny);
            log_info("Maximum source position coordinates: %5.2f (x) : %5.2f (y)\n", srec_maxx, srec_maxy);

            /* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
            if (((srec_maxx < 0.0) || (srec_maxy < 0.0)) || ((srec_minx < 0.0) || (srec_miny < 0.0))) {
                log_fatal("Coordinate of at least one source location is outside the global grid.\n");
            }

            if ((srec_maxx > gv->NX * gv->DH * gv->NPROCX) || (srec_maxy > gv->NY * gv->DH * gv->NPROCY)) {
                log_fatal("Coordinate of at least one source location is outside the global grid.\n");
            }

            /* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
            if ((srec_maxx < ((float)gv->FW * gv->DH)) || (srec_minx < ((float)gv->FW * gv->DH))) {
                /* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary") */
                log_warn
                    ("Coordinate of at least one source location is inside the Absorbing Boundary (left boundary).\n");
            }
            if (srec_maxx > (gv->NX * gv->DH * gv->NPROCX - (float)gv->FW * gv->DH)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one source location is inside the Absorbing Boundary (right boundary).\n");
            }
            if (srec_maxy > (gv->NY * gv->DH * gv->NPROCY - (float)gv->FW * gv->DH)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one source location is inside the Absorbing Boundary (lower boundary).\n");
            }
            if ((srec_miny < ((float)gv->FW * gv->DH)) && !(gv->FREE_SURF)) {
                /* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary") */
                log_warn
                    ("Coordinate of at least one source location is inside the Absorbing Boundary (top boundary).\n");
            }
        }

        log_info("------------------------- Absorbing frame checks ------------\n");
        if (gv->ABS_TYPE == 1) {
            log_info("CPML boundary (ABS_TYPE=1) with width of %d grid points (%5.2fm).\n", gv->FW,
                     (float)gv->FW * gv->DH);
            if (gv->FW < 10) {
                log_warn("Width (FW) of absorbing frame should be at least 10 grid points.\n");
                log_warn("Be aware of artificial reflections from grid boundaries!\n");
            }
        }
        if (gv->ABS_TYPE == 2) {
            log_info("Absorbing boundary (ABS_TYPE=2) with width of %d grid points (%5.2fm).\n", gv->FW,
                     (float)gv->FW * gv->DH);
            if (gv->FW < 30) {
                log_warn("Width (FW) of absorbing frame should be at least 30 grid points.\n");
                log_warn("Be aware of artificial reflections from grid boundaries!\n");
            }
        }
        if (((gv->NX) < gv->FW) || ((gv->NY) < gv->FW)) {
            log_error
                ("Width of absorbing boundary (FW=%d grid points) is larger than at least one subdomain dimension.\n",
                 gv->FW);
            log_error("Subdomain dimensions: NX/NPROCX=%d, NY/NPROCY=%d in grid points.\n", gv->NX, gv->NY);
            log_fatal
                ("You need to choose a smaller width of absorbing frame (FW) or increase the subdomain dimensions.\n");
        }
        if (gv->BOUNDARY) {
            if (gv->ABS_TYPE == 1 || gv->ABS_TYPE == 2) {
                log_warn("You have activated a periodic boundary and set an absorbing boundary at the same time!\n");
            }
        }
    }
}
