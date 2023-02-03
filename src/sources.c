
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

/* ----------------------------------------------------------------------
 * Reading (distributed) source positions, timeshift, centre frequency
 * and amplitude from SOURCE_FILE.
 * ---------------------------------------------------------------------- */

#include "fd.h"
#include "logging.h"

float **sources(int *nsrc, GlobVar * gv)
{
    float **srcpos = NULL;
    int i, l, isrc = 0, current_source = 0, nvarin = 0;
    float xsrc, ysrc, tshift, tan_phi, dz, x, fc = 0.0;
    char buffer[STRING_SIZE], bufferstring[10], cline[256];
    FILE *fpsrc = NULL;

    /* do NOT remove the FALLTHRU comments below, they are used to tell the compiler
     * (static code checker) that this is an intentional fall through */

    if (gv->MPID == 0) {
        if (gv->SRCREC == 1) {  /* read source positions from file */
            *nsrc = 0;
            log_info("------------------------- Source parameters (II) ------------\n");
            log_info("Reading source parameters from file %s.\n", gv->SOURCE_FILE);
            if ((fpsrc = fopen(gv->SOURCE_FILE, "r")) == NULL) {
                log_fatal("Source file %s could not be opened.\n", gv->SOURCE_FILE);
            }
            while (fgets(buffer, STRING_SIZE, fpsrc)) {
                sscanf(buffer, "%s", bufferstring);
                /* checks if the line contains a '#'character which indicates a comment line,
                 * and if the reading of a string was successful, which is not the case for an empty line */
                if ((strchr(buffer, '#') == 0) && (sscanf(buffer, "%s", bufferstring) == 1))
                    ++(*nsrc);
            }

            rewind(fpsrc);

            if ((nsrc) == 0)
                log_warn("Could not determine number of sources parameter sets; assuming %d.\n", (*nsrc = 0));
            else
                log_info("Number of source positions: %d\n", *nsrc);

            srcpos = matrix(1, NSPAR, 1, *nsrc);

            /* memory for source position definition (Ricker, Fuchs-Mueller, sin**3 & from_File) */
            if (gv->SOURCE_SHAPE <= 4) {
                /* srcpos[1][l] = x position
                 * srcpos[2][l] = depth position
                 * srcpos[3][l] = horizontal position (always zero in 2D)
                 * srcpos[4][l] = time delay (source start time)
                 * srcpos[5][l] = centre frequency
                 * srcpos[6][l] = amplitude
                 * srcpos[7][l] = azimuth [°] (optional)
                 * srcpos[8][l] = SOURCE_TYPE (optional)
                 */
                l = 0;
                while (fgets(cline, 255, fpsrc)) {
                    sscanf(cline, "%s", bufferstring);
                    if ((strchr(cline, '#') == 0) && (sscanf(cline, "%s", bufferstring) == 1)) {
                        ++l;
                        if (l > *nsrc)
                            log_fatal("sources.c: buffer not large enough to store all sources - programming error.\n");
                        nvarin =
                            sscanf(cline, "%f%f%f%f%f%f%f", &xsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l],
                                   &srcpos[7][l], &srcpos[8][l]);
                        switch (nvarin) {
                          case 0:
                              xsrc = 0.0;
                              /* FALLTHRU */
                          case 1:
                              ysrc = 0.0;
                              /* FALLTHRU */
                          case 2:
                              if (gv->MPID == 0)
                                  log_error("No time shift defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 3:
                              if (gv->MPID == 0)
                                  log_error("No frequency defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 4:
                              if (gv->MPID == 0)
                                  log_error("No amplitude defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 5:
                              srcpos[7][l] = 0.0;
                              /* FALLTHRU */
                          case 6:
                              srcpos[8][l] = gv->SOURCE_TYPE;
                        }
                        if ((srcpos[8][l] != 4) && (nvarin > 5)) {
                            current_source = (int)srcpos[8][l];
                            if (gv->MPID == 0)
                                log_warn("SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n",
                                         l, current_source);
                        }
                        srcpos[1][l] = xsrc;
                        srcpos[2][l] = ysrc;
                        srcpos[3][l] = 0.0;
                        srcpos[4][l] = tshift;
                        fc = srcpos[5][l];
                    }
                }
            }
            /* memory for source position definition (Berlage) */
            else if (5 == gv->SOURCE_SHAPE) {
                /* srcpos[1][l] = x position
                 * srcpos[2][l] = depth position
                 * srcpos[3][l] = horizontal position (always zero in 2D)
                 * srcpos[4][l] = time delay (source start time)
                 * srcpos[5][l] = centre frequency
                 * srcpos[6][l] = amplitude
                 * srcpos[7][l] = azimuth [°] (optional)
                 * srcpos[8][l] = SOURCE_TYPE (optional)
                 * srcpos[9][l] = time exponent (Berlage only)
                 * srcpos[10][l] = exponential decay factor (Berlage only)
                 * srcpos[11][l] = initial phase angle [°] (Berlage only)
                 */
                l = 0;
                while (fgets(cline, 255, fpsrc)) {
                    sscanf(cline, "%s", bufferstring);
                    if ((strchr(cline, '#') == 0) && (sscanf(cline, "%s", bufferstring) == 1)) {
                        ++l;
                        if (l > *nsrc)
                            log_fatal("sources.c: buffer not large enough to store all sources - programming error.\n");
                        nvarin =
                            sscanf(cline, "%f%f%f%f%f%f%f%f%f%f", &xsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l],
                                   &srcpos[9][l], &srcpos[10][l], &srcpos[11][l], &srcpos[7][l], &srcpos[8][l]);
                        switch (nvarin) {
                          case 0:
                              xsrc = 0.0;
                              /* FALLTHRU */
                          case 1:
                              ysrc = 0.0;
                              /* FALLTHRU */
                          case 2:
                              if (gv->MPID == 0)
                                  log_error("No time shift defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 3:
                              if (gv->MPID == 0)
                                  log_error("No frequency defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 4:
                              if (gv->MPID == 0)
                                  log_error("No amplitude defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 5:
                              if (gv->MPID == 0)
                                  log_error("No time exponent defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 6:
                              if (gv->MPID == 0)
                                  log_error("No exponential decay factor defined for source %i in %s!\n", l,
                                            gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 7:
                              if (gv->MPID == 0)
                                  log_error("No initial phase angle [°] defined for source %i in %s!\n", l,
                                            gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 8:
                              srcpos[7][l] = 0.0;
                              /* FALLTHRU */
                          case 9:
                              srcpos[8][l] = gv->SOURCE_TYPE;
                        }
                        if ((srcpos[8][l] != 4) && (nvarin > 8)) {
                            current_source = (int)srcpos[8][l];
                            if (gv->MPID == 0)
                                log_warn("SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n",
                                         l, current_source);
                        }
                        srcpos[1][l] = xsrc;
                        srcpos[2][l] = ysrc;
                        srcpos[3][l] = 0.0;
                        srcpos[4][l] = tshift;
                        fc = srcpos[5][l];
                    }
                }
            }
            /* memory for source position definition (Klauder) */
            else if (6 == gv->SOURCE_SHAPE) {
                /* srcpos[1][l] = x position
                 * srcpos[2][l] = depth position
                 * srcpos[3][l] = horizontal position (always zero in 2D)
                 * srcpos[4][l] = time delay (source start time)
                 * srcpos[5][l] = centre frequency
                 * srcpos[6][l] = amplitude
                 * srcpos[7][l] = azimuth [°] (optional)
                 * srcpos[8][l] = SOURCE_TYPE (optional)
                 * srcpos[9][l] = minimum frequency (Klauder only)
                 * srcpos[10][l] = maximum frequency (Klauder only)
                 * srcpos[11][l] = sweep length [s] (Klauder only)
                 * srcpos[12][l] = width of the Klauder wavelet in number of centre periods (Klauder only)
                 */
                l = 0;
                while (fgets(cline, 255, fpsrc)) {
                    sscanf(cline, "%s", bufferstring);
                    if ((strchr(cline, '#') == 0) && (sscanf(cline, "%s", bufferstring) == 1)) {
                        ++l;
                        if (l > *nsrc)
                            log_fatal("sources.c: buffer not large enough to store all sources - programming error.\n");
                        nvarin =
                            sscanf(cline, "%f%f%f%f%f%f%f%f%f%f", &xsrc, &ysrc, &tshift, &srcpos[9][l], &srcpos[10][l],
                                   &srcpos[6][l], &srcpos[11][l], &srcpos[12][l], &srcpos[7][l], &srcpos[8][l]);
                        switch (nvarin) {
                          case 0:
                              xsrc = 0.0;
                              /* FALLTHRU */
                          case 1:
                              ysrc = 0.0;
                              /* FALLTHRU */
                          case 2:
                              if (gv->MPID == 0)
                                  log_error("No time shift defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 3:
                              if (gv->MPID == 0)
                                  log_error("No minimum frequency defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 4:
                              if (gv->MPID == 0)
                                  log_error("No maximum frequency defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 5:
                              if (gv->MPID == 0)
                                  log_error("No amplitude defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 6:
                              if (gv->MPID == 0)
                                  log_error("No sweep length [s] defined for source %i in %s!\n", l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 7:
                              if (gv->MPID == 0)
                                  log_error
                                      ("No width of the wavelet (in number of center periods) defined for source %i in %s!\n",
                                       l, gv->SOURCE_FILE);
                              log_fatal("Missing parameter in SOURCE_FILE!\n");
                              /* FALLTHRU */
                          case 8:
                              srcpos[7][l] = 0.0;
                              /* FALLTHRU */
                          case 9:
                              srcpos[8][l] = gv->SOURCE_TYPE;
                        }
                        if ((srcpos[8][l] != 4) && (nvarin > 8)) {
                            current_source = (int)srcpos[8][l];
                            if (gv->MPID == 0)
                                log_warn("SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n",
                                         l, current_source);
                        }
                        srcpos[1][l] = xsrc;
                        srcpos[2][l] = ysrc;
                        srcpos[3][l] = 0.0;
                        srcpos[4][l] = tshift;
                        fc = (srcpos[9][l] + srcpos[10][l]) / 2;
                        srcpos[5][l] = fc;
                    }
                }
            }

            fclose(fpsrc);

            /* Compute maximum centre frequency */
            for (l = 1; l <= *nsrc; l++)
                if (srcpos[5][l] > fc)
                    fc = srcpos[5][l];
            log_info("Maximum frequency found: %6.2fHz\n", gv->SOURCE_FILE, fc);
            gv->TS = 1.0 / fc;

            /* outputs all sources per each subdomain / node */

            if (gv->MPID == 0) {
                if (gv->RUN_MULTIPLE_SHOTS)
                    log_info("All sources will be modelled individually (RUN_MULTIPLE_SHOTS=1).\n");
                else
                    log_info("All sources will be modelled simultaneously (RUN_MULTIPLE_SHOTS=0).\n");
            }

        }                       // end gv->SRCREC==1
        else if (gv->SRCREC == 2) {
            if (gv->PLANE_WAVE_DEPTH > 0) { /* plane wave excitation */
                if (gv->SOURCE_SHAPE > 4) {
                    log_fatal
                        ("Plane wave is only implemented for Ricker, Fuchs-Mueller, sin^3 or an external wavelet! Change parameter SOURCE_SHAPE!\n");
                }
                log_info("------------------------- Source parameters (II) ------------\n");
                log_info("Computing source nodes for plane wave excitation.\n");
                log_info("Depth: %5.2fm; incidence angle: %5.2fdeg\n", gv->PLANE_WAVE_DEPTH, gv->PLANE_WAVE_ANGLE);

                tan_phi = tan(gv->PLANE_WAVE_ANGLE * PI / 180.0);

                dz = (float)gv->NXG * gv->DH * tan_phi;
                log_info("Maximum depth of plane wave: %5.2fm\n", gv->PLANE_WAVE_DEPTH + dz);
                if ((gv->PLANE_WAVE_DEPTH + dz) <= gv->NYG * gv->DH) {
                    *nsrc = gv->NXG;
                    srcpos = matrix(1, 8, 1, *nsrc);
                    isrc = 0;
                    for (i = 1; i <= gv->NXG; i++) {
                        isrc++;
                        x = (float)i *gv->DH;
                        srcpos[1][isrc] = x;
                        srcpos[2][isrc] = gv->PLANE_WAVE_DEPTH + (tan_phi * x);
                        srcpos[3][isrc] = 0.0;
                        srcpos[4][isrc] = 0.0;
                        srcpos[5][isrc] = 1.0 / gv->TS;
                        srcpos[6][isrc] = 1.0;
                        srcpos[7][isrc] = 0.0;
                        srcpos[8][isrc] = gv->SOURCE_TYPE;
                    }
                } else
                    log_fatal("Maximum depth of plane wave exceeds model depth.\n");
            } else
                log_fatal("SRCREC parameter specifies PLANE_WAVE excitation, but PLANE_WAVE_DEPTH<=0!\n");
        } else
            log_fatal("SRCREC parameter is invalid (SRCREC!=1 or SRCREC!=2)! No source parameters specified!\n");
    }                           // end gv->MPID==0

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(nsrc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(gv->TS), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (gv->MPID != 0)
        srcpos = matrix(1, NSPAR, 1, *nsrc);
    MPI_Bcast(&srcpos[1][1], (*nsrc) * 12, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (gv->MPID == 0) {
        if (*nsrc > 50)
            log_warn("The following table is quite large (%d lines); only printing the first 50 entries!\n", *nsrc);
        if (4 >= gv->SOURCE_SHAPE) {
            log_info("  Shot         x         y   tshift       fc      amp azimuth type\n");
        } else if (5 == gv->SOURCE_SHAPE) {
            log_info("  Shot         x         y   tshift       fc      amp    n   alpha    phi0 azimuth type\n");
        } else if (6 == gv->SOURCE_SHAPE) {
            log_info
                ("  Shot         x         y   tshift     fmin     fmax       fc     amp   tsweep    twave azimuth type\n");
        }

        int maxprint = *nsrc;
        if (maxprint > 50)
            maxprint = 50;

        for (l = 1; l <= maxprint; l++) {
            if (4 >= gv->SOURCE_SHAPE) {
                log_info("%6d %9.2f %9.2f %8.4f %8.2f %8.2f %7.2f %4d\n",
                         l, srcpos[1][l], srcpos[2][l], srcpos[4][l], srcpos[5][l], srcpos[6][l], srcpos[7][l],
                         (int)srcpos[8][l]);
            } else if (5 == gv->SOURCE_SHAPE) {
                log_info("%6d %9.2f %9.2f %8.4f %8.2f %8.2f %4.1f %7.2f %7.2f %7.2f %4d\n",
                         l, srcpos[1][l], srcpos[2][l], srcpos[4][l], srcpos[5][l], srcpos[6][l], srcpos[9][l],
                         srcpos[10][l], srcpos[11][l], srcpos[7][l], (int)srcpos[8][l]);
            } else if (6 == gv->SOURCE_SHAPE) {
                log_info("%6d %9.2f %9.2f %8.4f %8.2f %8.2f %8.2f %7.2f %8.2f %8.2f %7.2f %4d\n",
                         l, srcpos[1][l], srcpos[2][l], srcpos[4][l], srcpos[9][l], srcpos[10][l], srcpos[5][l],
                         srcpos[6][l], srcpos[11][l], srcpos[12][l], srcpos[7][l], (int)srcpos[8][l]);
            }
        }
    }                           // end gv->MPID==0

    return srcpos;
}
