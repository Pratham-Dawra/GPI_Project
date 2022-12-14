
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
 *   Write FD-Parameters to stdout                           
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void write_par(GlobVar *gv)
{
    char file_ext[8];

    log_info("------------------------- Processors ------------------------\n");
    log_info("Number of subdomains in x-direction (NPROCX): %d\n", gv->NPROCX);
    log_info("Number of subdomains in y-direction (NPROCY): %d\n", gv->NPROCY);
    log_info("Total number of MPI processes/subdomains (NP): %d\n", gv->NP);

    log_info("------------------------- FD Algorithm ----------------------\n");
    log_info("Grid: standard staggered grid (SSG) (Virieux-grid)\n");
    log_info("Order of spatial FD operators (FDORDER): %d\n", gv->FDORDER);
    log_info("Order of temporal FD operator (FDORDER_TIME): %d\n", gv->FDORDER_TIME);

    log_info("------------------------- Wave Equation ---------------------\n");
    log_info("Type (WEQ): %s\n", get_weq_verbose(gv->WEQ));

    log_info("------------------------- Discretization --------------------\n");
    log_info("Model size (grid points) in x-direction (NX): %d\n", gv->NX);
    log_info("Model size (grid points) in y-direction (NY): %d\n", gv->NY);
    log_info("Grid-spacing (DH): %em\n", gv->DH);
    log_info("Model size (real) in x-direction: %em\n", (gv->NX - 1) * gv->DH);
    log_info("Model size (real) in y-direction: %em\n", (gv->NY - 1) * gv->DH);
    log_info("Time of wave propagation (TIME): %es\n", gv->TIME);
    log_info("Time step interval (DT): %es\n", gv->DT);
    log_info("Number of time steps: %d \n", gv->NT);

    log_info("------------------------- Source parameters -----------------\n");
    if (gv->SRCREC) {
        log_info("Source parameters will be read from ASCII file %s.\n", gv->SOURCE_FILE);
    } else {
        log_info("Plane wave excitation depth (PLANE_WAVE_DEPTH): %5.2fm\n", gv->PLANE_WAVE_DEPTH);
        log_info("Incidence angle of plane P-wave (from vertical; PLANE_WAVE_ANGLE): %5.2fdeg\n", gv->PLANE_WAVE_ANGLE);
        log_info("Duration of source signal (TS): %es\n", gv->TS);
        log_info("Center frequency is approximately %eHz.\n", 1.0 / gv->TS);
    }
    switch (gv->SOURCE_SHAPE) {
      case 1:
          log_info("Source wavelet: Ricker (zero-phase; shifted)\n");
          break;
      case 2:
          log_info("Source wavelet: Fuchs-Mueller\n");
          break;
      case 3:
          log_info("Source wavelet will be read from ASCII file %s.\n", gv->SIGNAL_FILE);
          break;
      case 4:
          log_info("Source wavelet: sine^3\n");
          break;
      case 5:
          log_info("Source wavelet: Berlage (minimum-phase)\n");
          break;
      case 6:
          log_info("Source wavelet: Klauder (zero-phase; shifted)\n");
          break;
      default:
          log_fatal("Sorry, incorrect specification of parameter SOURCE_SHAPE (source wavelet)!\n");
          break;
    }
    switch (gv->SOURCE_TYPE) {
      case 1:
          log_info("Type of source: explosive source\n");
          break;
      case 2:
          log_info("Type of source: point source with directive force in (horizontal) x-direction\n");
          break;
      case 3:
          log_info("Type of source: point source with directive force in (vertical) y-direction\n");
          break;
      case 4:
          log_info("Type of source: point source with directive force in (horizontal) z-direction\n");
          break;
      default:
          log_fatal("Sorry, incorrect specification of parameter SOURCE_TYPE (source type)!\n");
          break;
    }
    if (1 == gv->SIGOUT) {
        switch (gv->SIGOUT_FORMAT) {
          case 1:
              sprintf(file_ext, "su");
              break;
          case 2:
              sprintf(file_ext, "txt");
              break;
          case 3:
              sprintf(file_ext, "bin");
              break;
          default:
              log_fatal
                  ("Sorry, incorrect specification of parameter SIGOUT_FORMAT (output format of source signature)!\n");
        }
        log_info("Source signature will be written to file %s.%s.\n", gv->SIGOUT_FILE, file_ext);
    }

    if (gv->SEISMO) {
        log_info("------------------------- Receivers -------------------------\n");
        if (gv->READREC) {
            log_info("Receiver positions will be read from ASCII file %s.\n", gv->REC_FILE);
            log_info("Reference point for receivers (REFRECX, REFRECY): (%f, %f)m\n", gv->REFREC[1], gv->REFREC[2]);
        } else if (gv->REC_ARRAY > 0) {
            log_info("Using %d horizontal line(s) of receivers (REC_ARRAY).\n", gv->REC_ARRAY);
            log_info("Depth of upper line (REC_ARRAY_DEPTH): %em\n", gv->REC_ARRAY_DEPTH);
            log_info("Vertical increment between lines (REC_ARRAY_DIST): %em\n", gv->REC_ARRAY_DIST);
            log_info("Receiver distance (grid points) within line in x-direction (DRX): %d\n", gv->DRX);
        } else {
            log_info("First receiver position (XREC1, YREC1): (%e, %e)m\n", gv->XREC1, gv->YREC1);
            log_info("Last receiver position (XREC2, YREC2): (%e, %e)m\n", gv->XREC2, gv->YREC2);
        }
    }

    log_info("------------------------- Free surface ----------------------\n");
    if (gv->FREE_SURF)
        log_info("Free surface at the top of the model (FREE_SURF).\n");
    else
        log_info("No free surface at the top of the model (FREE_SURF).\n");

    log_info("------------------------- Absorbing frame -------------------\n");
    if (gv->FW > 0) {
        log_info("Width of absorbing frame (FW): %d grid points (%5.2fm)\n", gv->FW, (float)(gv->FW - 1) * gv->DH);
        if (gv->ABS_TYPE == 1) {
            log_info("Type of frame (ABS_TYPE): CPML\n");
            log_info("Damping velocity in the PML frame (VPPML): %fm/s\n", gv->VPPML);
            log_info("Frequency within the PML frame (FPML): %fHz\n", gv->FPML);
            log_info("NPOWER: %f\n", gv->NPOWER);
            log_info("K_MAX: %f\n", gv->K_MAX_CPML);
        }
        if (gv->ABS_TYPE == 2) {
            log_info("Type of frame (ABS_TYPE): exponential damping\n");
            log_info("Percentage of amplitude decay (DAMPING): %f\n", gv->DAMPING);
        }
    } else
        log_info("No absorbing frame used.\n");
    switch (gv->BOUNDARY) {
      case 0:
          log_info("No periodic boundary condition used (BOUNDARY).\n");
          break;
      case 1:
          log_info("Periodic boundary condition at left and right edges (BOUNDARY).\n");
          break;
      default:
          log_warn("Unknown value %d for parameter BOUNDARY; setting BOUNDARY=0 (no periodic boundary).", gv->BOUNDARY);
          gv->BOUNDARY = 0;
          break;
    }

    log_info("------------------------- Q approximation -------------------\n");
    log_info("Number of relaxation mechanisms (L): %d\n", gv->L);
    if (gv->L > 0) {
        log_info("Value for tau (TAU): %f\n", gv->TAU);
        for (int l = 1; l <= gv->L; l++) {
            log_info("Relaxation frequency %d (FL%d): %fHz\n", l, l, gv->FL[l]);
        }
    }

    if (gv->SNAP) {
        log_info("------------------------- Snapshots -------------------------\n");
        log_info("Snapshots: ");
        switch (gv->SNAP) {
          case 0:
              log_std("no shapshots will be output\n");
              break;
          case 1:
              log_std("vx, vy\n");
              break;
          case 2:
              log_std("p\n");
              break;
          case 3:
              log_std("curl, div\n");
              break;
          case 4:
              log_std("vx, vy, p, curl, div\n");
              break;
          default:
              log_fatal("Sorry, incorrect specification of parameter SNAP (snapshot output)!\n");
              break;
        }
        log_info("First snapshot (TSNAP1): %8.5fs\n", gv->TSNAP1);
        log_info("Last snapshot (TSNAP2): %8.5fs\n", gv->TSNAP2);
        log_info("Snapshot increment (TSNAPINC): %8.5fs\n", gv->TSNAPINC);
        log_info("First and last_horizontal grid point: %d, %d\n", 1, gv->NX);
        log_info("First and last vertical grid point: %d, %d\n", 1, gv->NY);
        log_info("Snapshot output file (SNAP_FILE): %s\n", gv->SNAP_FILE);
        switch (gv->SNAP_FORMAT) {
          case 1:
              log_fatal("SNAP_FORMAT=1, i.e. SU format, not yet available!\n");
              break;
          case 2:
              log_info("Output format (SNAP_FORMAT): ASCII\n");
              break;
          case 3:
              log_info("Output format (SNAP_FORMAT): binary (32-bit IEEE, native endian)\n");
              break;
          default:
              log_fatal("Sorry, incorrect specification of parameter SNAP_FORMAT (snapshot output format)!\n");
              break;
        }
    }

    if (gv->SEISMO) {
        log_info("------------------------- Seismograms -----------------------\n");
        switch (gv->SEIS_FORMAT) {
          case 1:
              sprintf(file_ext, "su");
              break;
          case 2:
              sprintf(file_ext, "txt");
              break;
          case 3:
              sprintf(file_ext, "bin");
              break;
          default:
              log_fatal("Sorry, incorrect specification of parameter SEIS_FORMAT (output format of seismograms)!\n");
              break;
        }
        log_info("Base name for output files (SEIS_FILE): %s\n", gv->SEIS_FILE);
        if ((gv->SEISMO == 1) || (gv->SEISMO == 4)) {
            log_info("Seismograms of vx and vy: %s_vx.%s, %s_vy.%s\n", gv->SEIS_FILE, file_ext, gv->SEIS_FILE,
                     file_ext);
        }
        if ((gv->SEISMO == 2) || (gv->SEISMO == 4)) {
            log_info("Seismograms of p: s_p.%s\n", gv->SEIS_FILE, file_ext);
        }
        if ((gv->SEISMO == 3) || (gv->SEISMO == 4)) {
            log_info("Seismograms of curl and div: t%s_rot.%s, %s_div.%s\n", gv->SEIS_FILE, file_ext, gv->SEIS_FILE,
                     file_ext);
        }
        switch (gv->SEIS_FORMAT) {
          case 1:
              log_info("Output format (SEIS_FORMAT): SU format (32-bit IEEE, native endian)\n");
              break;
          case 2:
              log_info("Output format (SEIS_FORMAT): ASCII\n");
              break;
          case 3:
              log_info("Output format (SEIS_FORMAT): binary (32-bit IEEE, native endian)\n");
              break;
          default:
              log_fatal("Sorry, incorrect specification of parameter SEIS_FORMAT (seismogram output format)!\n");
              break;
        }
        log_info("Sampling interval of output data: %fs\n", gv->NDT * gv->DT);
        if (!gv->READREC)
            log_info("Trace-spacing: %5.2fm\n", gv->NGEOPH * gv->DH);
        log_info("Number of samples per output trace: %d\n", iround((gv->NT - 1) / gv->NDT));
    }

    log_info("-------------------------------------------------------------\n");
    return;
}
