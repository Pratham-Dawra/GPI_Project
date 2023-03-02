
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

/*  ----------------------------------------------------------------------
 *  This is program SOFI2D.
 *  Parallel 2-D Viscoelastic Finite Difference Seismic Modeling *  using the Standard Staggered Grid (SSG)
 *
 *  PLEASE DO NOT DISTRIBUTE. PLEASE REFER OTHER PEOPLE TO :
 *
 *  Prof. Dr. Thomas Bohlen, Karlsruhe Institute of Technology,
 *  Geophysical Institute,
 *  Hertzstr. 16, 76187 Karlsruhe, Germany
 *  Phone/Fax: +49 (0)721 608 44416
 *  mailto:thomas.bohlen@kit.edu,
 *  http://www.gpi.kit.edu/
 *  http://www.gpi.kit.edu/SOFI2D.php
 *
 *  If you want to publish synthetic data calculated with this program please
 *  give a reference to the following paper:
 *  Bohlen, T., 2002, Parallel 3-D viscoelastic finite-difference seismic modelling,
 *  Computers @ Geopsciences, Vol. 28, No. 8, 887-889.
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar_struct.h"
#include "logging.h"
#include "macros.h"
#include "enums.h"
#ifdef EBUG
#include "debug_buffers.h"
#endif

#include <unistd.h>

int main(int argc, char **argv)
{
    /* variables in main */
    int nt;
    int lsnap, nsnap = 0;
    int nsrc_loc = 0;
//    int infocounter = 0;
    int ishot, nshots;          /* Added ishot and nshots for multiple shots */
    int **dummy = NULL;
    /*Limits for local grids defined in subgrid_bounds.c */
    char sigf[STRING_SIZE * 2], file_ext[5];
    clock_t cpu_time1 = 0, cpu_time = 0;
    FILE *log_fp = NULL;

    char ext[10];
    double time1 = 0.0, time2 = 0.0, time9 = 0.0;
//    double time_av_v_update = 0.0, time_av_s_update = 0.0,
//        time_av_v_exchange = 0.0, time_av_s_exchange = 0.0, time_av_timestep = 0.0;

    float *hc = NULL;
    float **signals = NULL;

    /* declare struct for global variables */
    GlobVar gv = {.MPID = -1,.OUTNTIMESTEPINFO = 100,.NDT = 1,.IDX = 1,.IDY = 1 };

    /* declare struct for acquisition variables */
    AcqVar acq = { };

    /* declare struct for wavefield variables */
    MemWavefield mpw = { };

    /* declare struct for model variables */
    MemModel mpm = { };

    /* declare struct for performance measures */
    Perform perf = { };

    /* initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &gv.NP);
    MPI_Comm_rank(MPI_COMM_WORLD, &gv.MPID);

    /* initialize logging */
    log_init(NULL);
    log_banner(LOG_SOFI);

    time1 = MPI_Wtime();
    cpu_time1 = clock();

    if (gv.MPID == 0) {
        if (argc != 2) {
            log_fatal
                ("Unexpected number of commandline arguments; single argument required: name of json parameter file.\n");
        }

        const char *fileinp = argv[1];

        /* check if parameter file can be opened */
        if (access(fileinp, R_OK) != 0) {
            log_fatal("Cannot open/read json parameter file %s.\n", fileinp);
        }

        /* check suffix of parameter file */
        if (!STRSTRCOMP(fileinp, ".json")) {
            log_fatal("Parameter file %s has no .json suffix.\n", fileinp);
        }

        /* read json parameter file */
        read_par_json(fileinp, &gv);
    }

    /* exchange parameters between MPI processes */
    exchange_par(&gv);

    /* set logging verbosity */
    log_set_level_from_string(&(gv.LOG_VERBOSITY[0]));

    /* check file system/output directories */
    check_fs(&gv);

    sprintf(ext, ".%i", gv.MPID);
    strcat(gv.LOG_FILE, ext);

    /* set up logging output */
    switch (gv.LOG) {
      case 0:                  /* logging to stdout/stderr for all ranks */
          log_fp = NULL;
          log_set_output(NULL);
          log_infoc(0, "Log messages sent to stdout/stderr on all MPI ranks.\n");
          break;
      case 1:                  /* logging to file for all ranks */
          log_infoc(0, "Now redirecting log messages to log file on all MPI ranks.\n");
          if ((log_fp = fopen(gv.LOG_FILE, "w")) == NULL) {
              log_fatal("Opening log file %s for writing failed.\n", gv.LOG_FILE);
          }
          log_set_output(log_fp);
          log_info("This is the log file %s generated by PE %d.\n", gv.LOG_FILE, gv.MPID);
          break;
      case 2:                  /* logging to stdout/stderr on rank 0, logging to file for all other ranks */
          if (0 == gv.MPID) {
              log_fp = NULL;
              log_set_output(NULL);
              log_info("Now redirecting log messages to log file on all MPI ranks except rank 0.\n");
          } else {
              if ((log_fp = fopen(gv.LOG_FILE, "w")) == NULL) {
                  log_fatal("Opening log file %s for writing failed.\n", gv.LOG_FILE);
              }
              log_set_output(log_fp);
              log_info("This is the log file %s generated by PE %d.\n", gv.LOG_FILE, gv.MPID);
          }
          break;
      default:
          log_warn("Unknown value %d for parameter LOG encountered; using LOG=0.\n", gv.LOG);
          gv.LOG = 0;
          log_fp = NULL;
          log_set_output(NULL);
          log_infoc(0, "Log messages sent to stdout/stderr on all MPI ranks.\n");
          break;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* domain decomposition */
    initproc(&gv);

    gv.NT = iround(gv.TIME / gv.DT);    /* number of time steps */
    gv.NS = iround(gv.NT / gv.NDT); /* number of samples per trace */
    lsnap = iround(gv.TSNAP1 / gv.DT);  /* first snapshot at this time step */

    /* output of parameters */
    if (gv.MPID == 0) {
        write_par(&gv);
    }

    /* Reading acquisition parameters */
    acq_read(&acq, &gv);

    /* memory allocation of buffers */
    initmem(&mpm, &mpw, &gv);

    /* initialize FD operators */
    initfd(&gv);

    /* Holberg coefficients for FD operators */
    hc = holbergcoeff(&gv);

    MPI_Barrier(MPI_COMM_WORLD);

    /* create model grids */
    readmod(&mpm, &gv);

    MPI_Barrier(MPI_COMM_WORLD);

    /* check if the FD run will be stable and free of numerical dispersion */
    switch (gv.WEQ) {
      case AC_ISO:             /* acoustic */
          log_fatal("not yet implemented\n");
          break;
      case AC_VTI:             /* acoustic VTI */
          log_fatal("not yet implemented\n");
          break;
      case AC_TTI:             /* acoustic TTI */
          log_fatal("not yet implemented\n");
          break;
      case EL_ISO:             /* elastic */
          checkfd(mpm.prho, mpm.ppi, mpm.pu, mpm.ptaus, mpm.ptaup, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos, &gv);
          break;
      case VEL_ISO:            /* viscoelastic */
          checkfd(mpm.prho, mpm.ppi, mpm.pu, mpm.ptaus, mpm.ptaup, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos, &gv);
          break;
      case EL_VTI:             /* elastic VTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptaus, mpm.ptaup, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos,
                  &gv);
          break;
      case VEL_VTI:            /* viscoelastic VTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptau55, mpm.ptau11, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos,
                  &gv);
          break;
      case EL_TTI:             /* elastic TTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptaus, mpm.ptaup, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos,
                  &gv);
          break;
      case VEL_TTI:            /* viscoelastic TTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptau55, mpm.ptau11, mpm.peta, hc, acq.srcpos, acq.nsrc, acq.recpos,
                  &gv);
          break;
      case VAC_ISO:            /* viscoacoustic */
          log_fatal("not yet implemented\n");
          break;
      case VAC_VTI:            /* viscoacoustic VTI */
          log_fatal("not yet implemented\n");
          break;
      case VAC_TTI:            /* viscoacoustic TTI */
          log_fatal("not yet implemented\n");
          break;
      default:
          log_fatal("Unknown WEQ.\n");
    }

    /* calculate damping coefficients for CPMLs */
    if (gv.ABS_TYPE == 1) {
        PML_pro(&mpm, &gv);
    }

    /* calculate 2-D array for exponential damping of reflections at the edges of the numerical mesh */
    if (gv.ABS_TYPE == 2) {
        absorb(mpm.absorb_coeff, &gv);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    prepmod(&mpm, &gv);

    MPI_Barrier(MPI_COMM_WORLD);

    time2 = MPI_Wtime();
    log_infoc(0, "Starting time stepping around real time %4.2fs.\n", time2 - time1);

    /*----------------------  loop over multiple shots  ------------------*/

    if (gv.RUN_MULTIPLE_SHOTS) {
        nshots = acq.nsrc;
    } else {
        nshots = 1;
    }

    for (ishot = 1; ishot <= nshots; ishot++) {

        for (nt = 1; nt <= 12; nt++) {
            acq.srcpos_current[nt][1] = acq.srcpos[nt][ishot];
        }

        if (gv.RUN_MULTIPLE_SHOTS) {
            log_info("Starting simulation for shot %d of %d.\n", ishot, nshots);
            //log_info("number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t    source_azimuth\n");
            //log_info("   %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f  \t %6.2f\n", ishot, acq.srcpos_current[1][1], acq.srcpos_current[2][1],
            //         acq.srcpos_current[4][1], acq.srcpos_current[5][1], acq.srcpos_current[6][1], acq.srcpos_current[7][1]);

            /* find this single source positions on subdomains  */
            if (nsrc_loc > 0)
                free_matrix(acq.srcpos_loc, 1, NSPAR, 1, 1);
            acq.srcpos_loc = splitsrc(acq.srcpos_current, &nsrc_loc, 1, &gv);
        } else {
            acq.srcpos_loc = splitsrc(acq.srcpos, &nsrc_loc, acq.nsrc, &gv);    /* Distribute source positions on subdomains */
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* calculate wavelet for each source point */
        signals = wavelet(acq.srcpos_loc, nsrc_loc, &gv);

        /* write source wavelet to file in subdomain that contains source */
        if (signals && 1 == gv.SIGOUT) {
            switch (gv.SIGOUT_FORMAT) {
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
                  log_fatal("Unknown SIGOUT_FORMAT encountered.\n");
                  break;
            }
            dummy = imatrix(1, 3, 1, 1);
            dummy[1][1] = iround(acq.srcpos_loc[1][ishot] / gv.DH);
            dummy[2][1] = iround(acq.srcpos_loc[2][ishot] / gv.DH);
            dummy[3][1] = 0;
            sprintf(sigf, "%s.shot%d.%s", gv.SIGOUT_FILE, ishot, file_ext);
            log_info("Writing source wavelet to file %s.\n", sigf);
            outseis_glob(fopen(sigf, "w"), signals, dummy, 1, acq.srcpos_loc, gv.NT, gv.SIGOUT_FORMAT, ishot, 0, &gv);
            free_imatrix(dummy, 1, 3, 1, 1);
        }

        /* initialize wavefield with zero */
        zero_wavefield(&mpw, &gv);

        subgrid_bounds(1, gv.NX, 1, gv.NY, &gv);

        /*---------------------------------------------------------------*/

        /*----------------------  loop over timesteps  ------------------*/

        /*---------------------------------------------------------------*/

        time_loop(ishot, lsnap, nsnap, hc, signals, nsrc_loc, &acq, &mpm, &mpw, &gv, &perf);

        /*---------------------------------------------------------------*/

        /*--------------------  End  of loop over timesteps ----------*/

        /*---------------------------------------------------------------*/

        saveseis(ishot, &acq, &gv);

    }   /*----------------------  end of loop over multiple shots  ------------------*/

    /* deallocate memory */
    freemem(&mpm, &mpw, &gv);

    if (gv.SEISMO)
        free_imatrix(acq.recpos, 1, 3, 1, gv.NTRG);

    MPI_Barrier(MPI_COMM_WORLD);

    if (gv.MPID == 0) {
        perf.time_av_v_update = perf.time_av_v_update / (double)perf.infocounter;
        perf.time_av_s_update = perf.time_av_s_update / (double)perf.infocounter;
        perf.time_av_v_exchange = perf.time_av_v_exchange / (double)perf.infocounter;
        perf.time_av_s_exchange = perf.time_av_s_exchange / (double)perf.infocounter;
        perf.time_av_timestep = perf.time_av_timestep / (double)perf.infocounter;
        log_info("Approximate average times for\n");
        log_info("  velocity update: .. %.6lfs.\n", perf.time_av_v_update);
        log_info("  stress update: .... %.6lfs.\n", perf.time_av_s_update);
        log_info("  velocity exchange:  %.6lfs.\n", perf.time_av_v_exchange);
        log_info("  stress exchange: .. %.6lfs.\n", perf.time_av_s_exchange);
        log_info("  time step: ........ %.6lfs.\n", perf.time_av_timestep);
        cpu_time = clock() - cpu_time1;
        log_info("CPU time of program per PE: %.3lfs.\n", (double)cpu_time / (double)CLOCKS_PER_SEC);
        time9 = MPI_Wtime();
        log_info("Total real time of program: %.3lfs.\n", time9 - time1);
    }

    /* finalize logging */
    log_finalize();

    if (log_fp)
        fclose(log_fp);

    /* finalize MPI */
    MPI_Finalize();

    return EXIT_SUCCESS;
}
