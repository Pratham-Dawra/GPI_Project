
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
 *
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
    int lsnap, nsnap = 0, lsamp = 0;
    int nsrc = 0, nsrc_loc = 0, infocounter = 0;
    int ishot, nshots;          /* Added ishot and nshots for multiple shots */
    int **dummy = NULL;
    /*Limits for local grids defined in subgrid_bounds.c */
    char sigf[STRING_SIZE * 2], file_ext[5];
    clock_t cpu_time1 = 0, cpu_time = 0;
    FILE *log_fp = NULL;

    char ext[10];
    double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0;
    double time5 = 0.0, time6 = 0.0, time7 = 0.0, time8 = 0.0;
    double time_av_v_update = 0.0, time_av_s_update = 0.0,
        time_av_v_exchange = 0.0, time_av_s_exchange = 0.0, time_av_timestep = 0.0;

    /* We need some pointers for the time shift for Adam Bashforth */
    float **shift_s1 = NULL, **shift_s2 = NULL;
    float **shift_v1 = NULL, **shift_v2 = NULL, **shift_v3 = NULL, **shift_v4 = NULL;
    float ***shift_r1 = NULL, ***shift_r2 = NULL, ***shift_r3 = NULL;

    float **srcpos = NULL, **srcpos_loc = NULL, **signals = NULL, *hc = NULL, **srcpos_current = NULL;
    int **recpos = NULL, **recpos_loc = NULL;

    int *recswitch = NULL;

    /* declare struct for global variables */
    GlobVar gv = {.MPID = -1,.OUTNTIMESTEPINFO = 100,.NDT = 1,.IDX = 1,.IDY = 1 };

    /* declare struct for wavefield variables */
    MemWavefield mpw = { };

    /* declare struct for model variables */
    MemModel mpm = { };

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
    lsamp = gv.NDT;

    /* output of parameters */
    if (gv.MPID == 0) {
        write_par(&gv);
    }

    /* NXG, NYG denote size of the entire (global) grid */
    gv.NXG = gv.NX;
    gv.NYG = gv.NY;

    /* In the following, NX and NY denote size of the local grid ! */
    gv.NX = gv.IENDX;
    gv.NY = gv.IENDY;

    /* memory allocation of buffers */
    initmem(&mpm, &mpw, &gv);

    /* initialize FD operators */
    initfd(&gv);

    /* Holberg coefficients for FD operators */
    hc = holbergcoeff(&gv);

    if (gv.SEISMO) {
        recpos = receiver(&gv);
        recswitch = ivector(1, gv.NTRG);
        recpos_loc = splitrec(recpos, recswitch, &gv);
    }

    /* memory for source position definition for saving the current positions */
    srcpos_current = matrix(1, NSPAR, 1, 1);

    /* Reading source positions from SOURCE_FILE */
    srcpos = sources(&nsrc, &gv);

    MPI_Barrier(MPI_COMM_WORLD);

    /* create model grids */
    if (gv.READMOD) {
        switch (gv.WEQ) {
          case AC_ISO:         /* acoustic */
              log_fatal("not yet implemented\n");
              break;
          case AC_VTI:         /* acoustic VTI */
              log_fatal("not yet implemented\n");
              break;
          case AC_TTI:         /* acoustic TTI */
              log_fatal("not yet implemented\n");
              break;
          case EL_ISO:         /* elastic */
              readmod_elastic(&mpm, &gv);
              break;
          case VEL_ISO:        /* viscoelastic */
              readmod_visco(&mpm, &gv);
              break;
          case EL_VTI:         /* elastic VTI */
              readmod_elastic_vti(&mpm, &gv);
#ifdef EBUG
              debug_check_matrix(mpm.prho, 0, gv.NX, gv.NY, 5, 0, "prho");
              debug_check_matrix(mpm.pc11, 0, gv.NX, gv.NY, 5, 0, "pc11");
              debug_check_matrix(mpm.pc33, 0, gv.NX, gv.NY, 5, 0, "pc33");
              debug_check_matrix(mpm.pc13, 0, gv.NX, gv.NY, 5, 0, "pc13");
              debug_check_matrix(mpm.pc55, 0, gv.NX, gv.NY, 5, 0, "pc55");
#endif
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              readmod_visco_vti(&mpm, &gv);
#ifdef EBUG
              debug_check_matrix(mpm.prho, 0, gv.NX, gv.NY, 6, 0, "prho");
              debug_check_matrix(mpm.pc11, 0, gv.NX, gv.NY, 6, 0, "pc11");
              debug_check_matrix(mpm.pc33, 0, gv.NX, gv.NY, 6, 0, "pc33");
              debug_check_matrix(mpm.pc13, 0, gv.NX, gv.NY, 6, 0, "pc13");
              debug_check_matrix(mpm.pc55, 0, gv.NX, gv.NY, 6, 0, "pc55");
              debug_check_matrix(mpm.ptau11, 0, gv.NX, gv.NY, 6, 0, "ptau11");
              debug_check_matrix(mpm.ptau33, 0, gv.NX, gv.NY, 6, 0, "ptau33");
              debug_check_matrix(mpm.ptau13, 0, gv.NX, gv.NY, 6, 0, "ptau13");
              debug_check_matrix(mpm.ptau55, 0, gv.NX, gv.NY, 6, 0, "ptau55");
              debug_check_vector(mpm.peta, 0, gv.L, 6, 0, "peta");
#endif
              break;
          case EL_TTI:         /* elastic TTI */
              readmod_elastic_tti(&mpm, &gv);
#ifdef EBUG
              debug_check_matrix(mpm.prho, 0, gv.NX, gv.NY, 7, 0, "prho");
              debug_check_matrix(mpm.pc11, 0, gv.NX, gv.NY, 7, 0, "pc11");
              debug_check_matrix(mpm.pc33, 0, gv.NX, gv.NY, 7, 0, "pc33");
              debug_check_matrix(mpm.pc13, 0, gv.NX, gv.NY, 7, 0, "pc13");
              debug_check_matrix(mpm.pc55, 0, gv.NX, gv.NY, 7, 0, "pc55");
              debug_check_matrix(mpm.pc15, 0, gv.NX, gv.NY, 7, 0, "pc15");
              debug_check_matrix(mpm.pc35, 0, gv.NX, gv.NY, 7, 0, "pc35");
#endif
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              readmod_visco_tti(&mpm, &gv);
#ifdef EBUG
              debug_check_matrix(mpm.prho, 0, gv.NX, gv.NY, 8, 0, "prho");
              debug_check_matrix(mpm.pc11, 0, gv.NX, gv.NY, 8, 0, "pc11");
              debug_check_matrix(mpm.pc33, 0, gv.NX, gv.NY, 8, 0, "pc33");
              debug_check_matrix(mpm.pc13, 0, gv.NX, gv.NY, 8, 0, "pc13");
              debug_check_matrix(mpm.pc55, 0, gv.NX, gv.NY, 8, 0, "pc55");
              debug_check_matrix(mpm.pc15, 0, gv.NX, gv.NY, 8, 0, "pc15");
              debug_check_matrix(mpm.pc35, 0, gv.NX, gv.NY, 8, 0, "pc35");
              debug_check_matrix(mpm.ptau11, 0, gv.NX, gv.NY, 8, 0, "ptau11");
              debug_check_matrix(mpm.ptau33, 0, gv.NX, gv.NY, 8, 0, "ptau33");
              debug_check_matrix(mpm.ptau13, 0, gv.NX, gv.NY, 8, 0, "ptau13");
              debug_check_matrix(mpm.ptau55, 0, gv.NX, gv.NY, 8, 0, "ptau55");
              debug_check_matrix(mpm.ptau15, 0, gv.NX, gv.NY, 8, 0, "ptau15");
              debug_check_matrix(mpm.ptau35, 0, gv.NX, gv.NY, 8, 0, "ptau35");
              debug_check_vector(mpm.peta, 0, gv.L, 8, 0, "peta");
#endif
              break;
          case VAC_ISO:        /* viscoacoustic */
              log_fatal("not yet implemented\n");
              break;
          case VAC_VTI:        /* viscoacoustic VTI */
              log_fatal("not yet implemented\n");
              break;
          case VAC_TTI:        /* viscoacoustic TTI */
              log_fatal("not yet implemented\n");
              break;
          default:
              log_fatal("Unknown WEQ.\n");
        }
    } else {
        switch (gv.WEQ) {
          case EL_ISO:         /* elastic */
              model_elastic(&mpm, &gv);
              break;
          case VEL_ISO:        /* viscoelastic */
              model_visco(&mpm, &gv);
              break;
          case EL_VTI:         /* elastic VTI */
              model_elastic_VTI(&mpm, &gv);
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              model_visco_vti(&mpm, &gv);
              break;
          case EL_TTI:         /* elastic TTI */
              model_elastic_TTI(&mpm, &gv);
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              model_visco_tti(&mpm, &gv);
              break;
          default:
              log_fatal("Internal model for your chosen WEQ not implemented.\n");
              break;
        }
    }

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
          checkfd(mpm.prho, mpm.ppi, mpm.pu, mpm.ptaus, mpm.ptaup, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
          break;
      case VEL_ISO:            /* viscoelastic */
          checkfd(mpm.prho, mpm.ppi, mpm.pu, mpm.ptaus, mpm.ptaup, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
          break;
      case EL_VTI:             /* elastic VTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptaus, mpm.ptaup, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
          break;
      case VEL_VTI:            /* viscoelastic VTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptau55, mpm.ptau11, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
          break;
      case EL_TTI:             /* elastic TTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptaus, mpm.ptaup, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
          break;
      case VEL_TTI:            /* viscoelastic TTI */
          checkfd(mpm.prho, mpm.pc11, mpm.pc55, mpm.ptau55, mpm.ptau11, mpm.peta, hc, srcpos, nsrc, recpos, &gv);
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

    /* For the calculation of the material parameters beteween gridpoints
     * they have to be averaged. For this, values lying at 0 and NX+1,
     * for example, are required on the local grid. These are now copied from the
     * neighbouring grids */

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
          matcopy_elastic(mpm.prho, mpm.ppi, mpm.pu, &gv);
          av_mue(mpm.pu, mpm.puipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
          break;
      case VEL_ISO:            /* viscoelastic */
          matcopy(mpm.prho, mpm.ppi, mpm.pu, mpm.ptaus, mpm.ptaup, &gv);    //???? why here just matcopy and not matcopy_elastic
          av_mue(mpm.pu, mpm.puipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
          av_tau(mpm.ptaus, mpm.ptausipjp, &gv);
          break;
      case EL_VTI:             /* elastic VTI */
          matcopy_elastic(mpm.prho, mpm.pc11, mpm.pc55, &gv);
          av_mue(mpm.pc55, mpm.pc55ipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
#ifdef EBUG
          debug_check_matrix(mpm.pc55ipjp, 0, gv.NX, gv.NY, 55, 0, "pc55ipjp");
          debug_check_matrix(mpm.prip, 0, gv.NX, gv.NY, 55, 0, "prip");
          debug_check_matrix(mpm.prjp, 0, gv.NX, gv.NY, 55, 0, "prjp");
#endif
          break;
      case VEL_VTI:            /* viscoelastic VTI */
          matcopy_elastic(mpm.prho, mpm.ptau55, mpm.pc55, &gv);
          av_mue(mpm.pc55, mpm.pc55ipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
          av_tau(mpm.ptau55, mpm.ptau55ipjp, &gv);
#ifdef EBUG
          debug_check_matrix(mpm.pc55ipjp, 0, gv.NX, gv.NY, 66, 0, "pc55ipjp");
          debug_check_matrix(mpm.prip, 0, gv.NX, gv.NY, 66, 0, "prip");
          debug_check_matrix(mpm.prjp, 0, gv.NX, gv.NY, 66, 0, "prjp");
          debug_check_matrix(mpm.ptau55ipjp, 0, gv.NX, gv.NY, 66, 0, "ptau55ipjp");
#endif
          break;
      case EL_TTI:             /* elastic TTI */
          matcopy_elastic(mpm.prho, mpm.pc11, mpm.pc55, &gv);
          matcopy_elastic(mpm.prho, mpm.pc15, mpm.pc35, &gv);
          av_mue(mpm.pc55, mpm.pc55ipjp, &gv);
          av_mue(mpm.pc15, mpm.pc15ipjp, &gv);
          av_mue(mpm.pc35, mpm.pc35ipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
#ifdef EBUG
          debug_check_matrix(mpm.pc55ipjp, 0, gv.NX, gv.NY, 77, 0, "pc55ipjp");
          debug_check_matrix(mpm.prip, 0, gv.NX, gv.NY, 77, 0, "prip");
          debug_check_matrix(mpm.prjp, 0, gv.NX, gv.NY, 77, 0, "prjp");
          debug_check_matrix(mpm.pc15ipjp, 0, gv.NX, gv.NY, 77, 0, "pc15ipjp");
          debug_check_matrix(mpm.pc35ipjp, 0, gv.NX, gv.NY, 77, 0, "pc35ipjp");
#endif
          break;
      case VEL_TTI:            /* viscoelastic TTI */
          matcopy_elastic(mpm.prho, mpm.ptau55, mpm.pc55, &gv);
          matcopy_elastic(mpm.prho, mpm.ptau15, mpm.pc15, &gv);
          matcopy_elastic(mpm.prho, mpm.ptau35, mpm.pc35, &gv);
          av_mue(mpm.pc55, mpm.pc55ipjp, &gv);
          av_mue(mpm.pc15, mpm.pc15ipjp, &gv);
          av_mue(mpm.pc35, mpm.pc35ipjp, &gv);
          av_rho(mpm.prho, mpm.prip, mpm.prjp, &gv);
          av_tau(mpm.ptau55, mpm.ptau55ipjp, &gv);
          av_tau(mpm.ptau15, mpm.ptau15ipjp, &gv);
          av_tau(mpm.ptau35, mpm.ptau35ipjp, &gv);
#ifdef EBUG
          debug_check_matrix(mpm.pc55ipjp, 0, gv.NX, gv.NY, 88, 0, "pc55ipjp");
          debug_check_matrix(mpm.prip, 0, gv.NX, gv.NY, 88, 0, "prip");
          debug_check_matrix(mpm.prjp, 0, gv.NX, gv.NY, 88, 0, "prjp");
          debug_check_matrix(mpm.pc15ipjp, 0, gv.NX, gv.NY, 88, 0, "pc15ipjp");
          debug_check_matrix(mpm.pc35ipjp, 0, gv.NX, gv.NY, 88, 0, "pc35ipjp");
          debug_check_matrix(mpm.ptau55ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau55ipjp");
          debug_check_matrix(mpm.ptau15ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau15ipjp");
          debug_check_matrix(mpm.ptau35ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau35ipjp");
#endif
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

    /* Preparing memory variables for update_s (viscoelastic only) */

    if (gv.FDORDER_TIME == 2) {
        switch (gv.WEQ) {
          case AC_ISO:         /* acoustic */
              break;
          case AC_VTI:         /* acoustic VTI */
              break;
          case AC_TTI:         /* acoustic TTI */
              break;
          case EL_ISO:         /* elastic */
              break;
          case VEL_ISO:        /* viscoelastic */
              prepare_update_s(&mpm, &gv);
              break;
          case EL_VTI:
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc11);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc33);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc13);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc55ipjp);
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              prepare_update_s_vti(&mpm, &gv);

              break;
          case EL_TTI:
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc11);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc33);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc13);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc15);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc35);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc55ipjp);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc35ipjp);
              dt_mult(gv.NX, gv.NY, gv.DT, mpm.pc15ipjp);
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              prepare_update_s_tti(&mpm, &gv);
              break;
          case VAC_ISO:        /* viscoacoustic */
              log_fatal("not yet implemented\n");
              break;
          case VAC_VTI:        /* viscoacoustic VTI */
              log_fatal("not yet implemented\n");
              break;
          case VAC_TTI:        /* viscoacoustic TTI */
              log_fatal("not yet implemented\n");
              break;
          default:
              log_fatal("Unknown WEQ.\n");
        }
    }

    if ((gv.WEQ == VEL_ISO) && gv.FDORDER_TIME == 4) {
        prepare_update_s_4(&mpm, &gv);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    time2 = MPI_Wtime();
    log_infoc(0, "Starting time stepping around real time %4.2fs.\n", time2 - time1);

    /*----------------------  loop over multiple shots  ------------------*/

    if (gv.RUN_MULTIPLE_SHOTS) {
        nshots = nsrc;
    } else {
        nshots = 1;
    }

    for (ishot = 1; ishot <= nshots; ishot++) {

        for (nt = 1; nt <= 12; nt++) {
            srcpos_current[nt][1] = srcpos[nt][ishot];
        }

        if (gv.RUN_MULTIPLE_SHOTS) {
            log_info("Starting simulation for shot %d of %d.\n", ishot, nshots);
            //log_info("number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t    source_azimuth\n");
            //log_info("   %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f  \t %6.2f\n", ishot, srcpos_current[1][1], srcpos_current[2][1],
            //         srcpos_current[4][1], srcpos_current[5][1], srcpos_current[6][1], srcpos_current[7][1]);

            /* find this single source positions on subdomains  */
            if (nsrc_loc > 0)
                free_matrix(srcpos_loc, 1, NSPAR, 1, 1);
            srcpos_loc = splitsrc(srcpos_current, &nsrc_loc, 1, &gv);
        } else {
            srcpos_loc = splitsrc(srcpos, &nsrc_loc, nsrc, &gv);    /* Distribute source positions on subdomains */
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* calculate wavelet for each source point */
        signals = wavelet(srcpos_loc, nsrc_loc, &gv);

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
            dummy[1][1] = iround(srcpos_loc[1][ishot] / gv.DH);
            dummy[2][1] = iround(srcpos_loc[2][ishot] / gv.DH);
            dummy[3][1] = 0;
            sprintf(sigf, "%s.shot%d.%s", gv.SIGOUT_FILE, ishot, file_ext);
            log_info("Writing source wavelet to file %s.\n", sigf);
            outseis_glob(fopen(sigf, "w"), signals, dummy, 1, srcpos_loc, gv.NT, gv.SIGOUT_FORMAT, ishot, 0, &gv);
            free_imatrix(dummy, 1, 3, 1, 1);
        }

        /* initialize wavefield with zero */

        if (gv.ABS_TYPE == 1) {
            if (gv.L)
                zero_PML_visc(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw, &gv);
            else
                zero_PML_elastic(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw, &gv);
        }

        if (gv.FDORDER_TIME == 4) {
            if (gv.L) {
                zero_visco_4(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw, &gv);
            } else {
                zero_elastic_4(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw);
            }
        }

        if (gv.ABS_TYPE != 1) {
            if (gv.L)
                zero_visc(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw, &gv);
            else
                zero_elastic(-gv.ND + 1, gv.NX + gv.ND, -gv.ND + 1, gv.NY + gv.ND, &mpw);
        }

        /* Reseting lsmap to NDT for saving seismograms  */
        lsamp = gv.NDT;

        subgrid_bounds(1, gv.NX, 1, gv.NY, &gv);

        /*---------------------------------------------------------------*/

        /*----------------------  loop over timesteps  ------------------*/

        /*---------------------------------------------------------------*/

        for (nt = 1; nt <= gv.NT; nt++) {

            if (isnan(mpw.pvy[gv.NY / 2][gv.NX / 2])) {
                log_error("Time step: %d; pvy: %f.\n", nt, mpw.pvy[gv.NY / 2][gv.NX / 2]);
                log_fatal("Simulation is unstable!\n");
            }

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                log_info("Computing time step %d of %d.\n", nt, gv.NT);
                time3 = MPI_Wtime();
            }

            /*---------------------------------------------------------------*/
            /* update of particle velocities -------------------------------- */

            /*---------------------------------------------------------------*/
            if (gv.FDORDER_TIME == 2) {
                update_v_interior(nt, srcpos_loc, signals, nsrc_loc, hc, &mpm, &mpw, &gv);
#ifdef EBUG
                debug_check_matrix(mpw.pvx, nt, gv.NX, gv.NY, 121, 0, "pvx");
                debug_check_matrix(mpw.pvy, nt, gv.NX, gv.NY, 121, 0, "pvy");
#endif

                if (gv.FW) {
                    if (gv.ABS_TYPE == 1) {
                        update_v_PML(gv.NX, gv.NY, nt, &mpm, &mpw, &gv);
                    }

                    if (gv.ABS_TYPE == 2) {
                        update_v_abs(&mpm, &mpw, &gv);
                    }
#ifdef EBUG
                    debug_check_matrix(mpw.pvx, nt, gv.NX, gv.NY, 122, 0, "pvx");
                    debug_check_matrix(mpw.pvy, nt, gv.NX, gv.NY, 122, 0, "pvy");
#endif
                }
            }

            if (gv.FDORDER_TIME == 4) {
                update_v_interior_4(nt, srcpos_loc, signals, nsrc_loc, hc, &mpm, &mpw, &gv);
                if (gv.FW) {
                    if (gv.ABS_TYPE == 1) {
                        update_v_PML_4(gv.NX, gv.NY, nt, &mpm, &mpw, &gv);
                    }

                    if (gv.ABS_TYPE == 2) {
                        update_v_abs_4(nt, &mpm, &mpw, &gv);
                    }
                }

                /* Shift spatial derivations of the stress one time step back */
                shift_s1 = mpw.svx_4;
                mpw.svx_4 = mpw.svx_3;
                mpw.svx_3 = mpw.svx_2;
                mpw.svx_2 = mpw.svx_1;
                mpw.svx_1 = shift_s1;
                shift_s2 = mpw.svy_4;
                mpw.svy_4 = mpw.svy_3;
                mpw.svy_3 = mpw.svy_2;
                mpw.svy_2 = mpw.svy_1;
                mpw.svy_1 = shift_s2;
            }

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                time4 = MPI_Wtime();
                time_av_v_update += (time4 - time3);
                log_debug("Starting particle velocity exchange between PEs...\n");
            }

            /*---------------------------------------------------------------*/
            /* ------- exchange of particle velocities between PEs -------------- */

            /*---------------------------------------------------------------*/

            exchange_v(mpw.pvx, mpw.pvy, &mpw, &gv);

            if ((gv.WEQ == EL_TTI) || (gv.WEQ == VEL_TTI)) {    /* TTI */
                v_derivatives(&mpw, &gv);
                exchange_v(mpw.pvxx, mpw.pvyy, &mpw, &gv);
                exchange_v(mpw.pvyx, mpw.pvxy, &mpw, &gv);
            }

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                time5 = MPI_Wtime();
                time_av_v_exchange += (time5 - time4);
                log_debug("Finished particle velocity exchange between PEs (real time: %.4fs).\n", time5 - time4);
            }

            /*---------------------------------------------------------------*/
            /* stress update ------------------------------------------------ */

            /*---------------------------------------------------------------*/
            if (gv.FDORDER_TIME == 2) {

                switch (gv.WEQ) {
                  case AC_ISO: /* acoustic */
                      log_fatal("not yet implemented\n");
                      break;
                  case AC_VTI: /* acoustic VTI */
                      log_fatal("not yet implemented\n");
                      break;
                  case AC_TTI: /* acoustic TTI */
                      log_fatal("not yet implemented\n");
                      break;
                  case EL_ISO: /* elastic */
                      update_s_elastic_interior(nt, hc, &mpm, &mpw, &gv);

                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_elastic_PML(nt, &mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_elastic_abs(nt, &mpm, &mpw, &gv);
                      }
                      break;
                  case VEL_ISO:    /* viscoelastic */
                      update_s_visc_interior(nt, &mpm, &mpw, &gv);
                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_visc_PML(nt, &mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_visc_abs(&mpm, &mpw, &gv);
                      }
                      break;
                  case EL_VTI: /* elastic VTI */
                      update_s_elastic_VTI_interior(nt, &mpm, &mpw, &gv);
                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_elastic_VTI_PML(&mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_elastic_vti_abs(&mpm, &mpw, &gv);
                      }
#ifdef EBUG
                      debug_check_matrix(mpw.psxx, nt, gv.NX, gv.NY, 555, 0, "psxx");
                      debug_check_matrix(mpw.psyy, nt, gv.NX, gv.NY, 555, 0, "psyy");
                      debug_check_matrix(mpw.psxy, nt, gv.NX, gv.NY, 555, 0, "psxy");
#endif
                      break;
                  case VEL_VTI:    /* viscoelastic VTI */
                      update_s_visc_VTI_interior(nt, &mpm, &mpw, &gv);
                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_visc_VTI_PML(&mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_visc_vti_abs(&mpm, &mpw, &gv);
                      }
#ifdef EBUG
                      debug_check_matrix(mpw.psxx, nt, gv.NX, gv.NY, 666, 0, "psxx");
                      debug_check_matrix(mpw.psyy, nt, gv.NX, gv.NY, 666, 0, "psyy");
                      debug_check_matrix(mpw.psxy, nt, gv.NX, gv.NY, 666, 0, "psxy");
#endif
                      break;
                  case EL_TTI: /* elastic TTI */
                      update_s_elastic_TTI_interior(nt, &mpm, &mpw, &gv);

                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_elastic_TTI_PML(&mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_elastic_tti_abs(&mpm, &mpw, &gv);
                      }
#ifdef EBUG
                      debug_check_matrix(mpw.psxx, nt, gv.NX, gv.NY, 777, 0, "psxx");
                      debug_check_matrix(mpw.psyy, nt, gv.NX, gv.NY, 777, 0, "psyy");
                      debug_check_matrix(mpw.psxy, nt, gv.NX, gv.NY, 777, 0, "psxy");
#endif
                      break;
                  case VEL_TTI:    /* viscoelastic TTI */
                      update_s_visc_TTI_interior(nt, &mpm, &mpw, &gv);
                      if (gv.FW) {
                          if (gv.ABS_TYPE == 1)
                              update_s_visc_TTI_PML(&mpm, &mpw, &gv);
                          if (gv.ABS_TYPE == 2)
                              update_s_visc_tti_abs(&mpm, &mpw, &gv);
                      }
#ifdef EBUG
                      debug_check_matrix(mpw.psxx, nt, gv.NX, gv.NY, 888, 0, "psxx");
                      debug_check_matrix(mpw.psyy, nt, gv.NX, gv.NY, 888, 0, "psyy");
                      debug_check_matrix(mpw.psxy, nt, gv.NX, gv.NY, 888, 0, "psxy");
#endif
                      break;
                  case VAC_ISO:    /* viscoacoustic */
                      log_fatal("not yet implemented\n");
                      break;
                  case VAC_VTI:    /* viscoacoustic VTI */
                      log_fatal("not yet implemented\n");
                      break;
                  case VAC_TTI:    /* viscoacoustic TTI */
                      log_fatal("not yet implemented\n");
                      break;
                  default:
                      log_fatal("Unknown WEQ.\n");
                }
            }

            if (gv.FDORDER_TIME == 4) {
                if (gv.L) {     /* viscoelastic */
                    /* Not supported right now */
                    update_s_visc_interior_4(nt, hc, &mpm, &mpw, &gv);
                    if (gv.FW) {
                        if (gv.ABS_TYPE == 1) {
                            update_s_visc_PML_4(nt, &mpm, &mpw, &gv);
                        }
                        if (gv.ABS_TYPE != 1) {
                            update_s_visc_abs_4(nt, &mpm, &mpw, &gv);
                        }
                    }
                    /* Shift memory variables one time step back */
                    shift_r1 = mpw.pp_4;
                    mpw.pp_4 = mpw.pp_3;
                    mpw.pp_3 = mpw.pp_2;
                    mpw.pp_2 = mpw.pp;
                    mpw.pp = shift_r1;
                    shift_r2 = mpw.pr_4;
                    mpw.pr_4 = mpw.pr_3;
                    mpw.pr_3 = mpw.pr_2;
                    mpw.pr_2 = mpw.pr;
                    mpw.pr = shift_r2;
                    shift_r3 = mpw.pq_4;
                    mpw.pq_4 = mpw.pq_3;
                    mpw.pq_3 = mpw.pq_2;
                    mpw.pq_2 = mpw.pq;
                    mpw.pq = shift_r3;
                } else {        /* elastic */
                    update_s_elastic_interior_4(nt, hc, &mpm, &mpw, &gv);

                    if (gv.FW) {
                        if (gv.ABS_TYPE == 1)
                            update_s_elastic_PML_4(nt, &mpm, &mpw, &gv);
                        if (gv.ABS_TYPE != 1)
                            update_s_elastic_abs_4(nt, &mpm, &mpw, &gv);
                    }
                }
                /* Shift spatial derivatives from the velocity one time step back */
                shift_v1 = mpw.vxx_4;
                mpw.vxx_4 = mpw.vxx_3;
                mpw.vxx_3 = mpw.vxx_2;
                mpw.vxx_2 = mpw.vxx_1;
                mpw.vxx_1 = shift_v1;
                shift_v2 = mpw.vyy_4;
                mpw.vyy_4 = mpw.vyy_3;
                mpw.vyy_3 = mpw.vyy_2;
                mpw.vyy_2 = mpw.vyy_1;
                mpw.vyy_1 = shift_v2;
                shift_v3 = mpw.vxy_4;
                mpw.vxy_4 = mpw.vxy_3;
                mpw.vxy_3 = mpw.vxy_2;
                mpw.vxy_2 = mpw.vxy_1;
                mpw.vxy_1 = shift_v3;
                shift_v4 = mpw.vyx_4;
                mpw.vyx_4 = mpw.vyx_3;
                mpw.vyx_3 = mpw.vyx_2;
                mpw.vyx_2 = mpw.vyx_1;
                mpw.vyx_1 = shift_v4;
            }

            /* explosive source */
            if (gv.SOURCE_TYPE == 1)
                psource(nt, srcpos_loc, signals, nsrc_loc, &mpw, &gv);

            if ((gv.FREE_SURF) && (gv.POS[2] == 0)) {
                if (gv.L)       /* viscoelastic */
                    surface(1, hc, &mpm, &mpw, &gv);
                else
                    /* elastic */
                    surface_elastic(1, hc, &mpm, &mpw, &gv);
            }

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                time6 = MPI_Wtime();
                time_av_s_update += (time6 - time5);
                log_debug("Starting stress exchange between PEs...\n");
            }

            /*---------------------------------------------------------------*/
            /* -------- stress exchange between PEs -------- */

            /*---------------------------------------------------------------*/
            exchange_s(&mpw, &gv);

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                time7 = MPI_Wtime();
                time_av_s_exchange += (time7 - time6);
                log_debug("Finished stress exchange between PEs (real time: %.4fs).\n", time7 - time6);
            }

            /* store amplitudes at receivers in section-arrays */
            if ((gv.SEISMO) && (nt == lsamp) && (nt < gv.NT)) {

                seismo_ssg(lsamp, recpos_loc, hc, &mpm, &mpw, &gv);
                lsamp += gv.NDT;
            }

            /* WRITE SNAPSHOTS TO DISK */
            if ((gv.SNAP) && (nt == lsnap) && (nt <= gv.TSNAP2 / gv.DT)) {

                snap(nt, ++nsnap, hc, &mpm, &mpw, &gv);
                lsnap = lsnap + iround(gv.TSNAPINC / gv.DT);
            }

            if ((gv.MPID == 0) && ((nt - 1) % gv.OUTNTIMESTEPINFO == 0)) {
                ++infocounter;
                time8 = MPI_Wtime();
                time_av_timestep += (time8 - time3);

                // when we reach this point, we have completed nt out of gv.NT time steps; however, the
                // time_av_timestep variable has only been updated every gv.OUTNTIMESTEPINFO time step;
                // use infocounter to calculate correct average
                log_info("Total real time for time step %d: %.4fs. Shot %d, time left: %.2lfs.\n", nt, time8 - time3,
                         ishot, (gv.NT - nt) * time_av_timestep / (double)infocounter);
            }
        }

        /*---------------------------------------------------------------*/

        /*--------------------  End  of loop over timesteps ----------*/

        /*---------------------------------------------------------------*/

        log_infoc(0, "Finished time stepping.\n");

        /* write seismograms to file(s) */
        if (gv.SEISMO) {

            /* saves seismograms portion of each PE individually to file */
            //if (gv.NTR> 0) saveseis(SECTIONVX,SECTIONVY,SECTIONP,SECTIONCURL,SECTIONDIV,recpos,recpos_loc,gv.NTR,srcpos_current,ishot,gv.NS);

            /* merge of seismogram data from all PE and output data collectively */
            switch (gv.SEISMO) {
              case 1:          /* particle velocities only */
                  catseis(gv.SECTIONVX, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 1, &gv);
                  catseis(gv.SECTIONVY, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 2, &gv);

                  break;
              case 2:          /* pressure only */
                  catseis(gv.SECTIONP, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 4, &gv);

                  break;
              case 3:          /* curl and div only */
                  catseis(gv.SECTIONDIV, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 5, &gv);
                  catseis(gv.SECTIONCURL, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 6, &gv);

                  break;
              case 4:          /* everything */
                  catseis(gv.SECTIONVX, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 1, &gv);
                  catseis(gv.SECTIONVY, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 2, &gv);
                  catseis(gv.SECTIONP, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 4, &gv);
                  catseis(gv.SECTIONDIV, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 5, &gv);
                  catseis(gv.SECTIONCURL, gv.SEISMO_FULLDATA, recswitch, gv.NTRG, gv.NS);
                  if (gv.MPID == 0)
                      saveseis_glob(gv.SEISMO_FULLDATA, recpos, srcpos, ishot, gv.NS, 6, &gv);

                  break;
              default:
                  break;

            }
        }
    }                           /* end of loop over shots */

    /* deallocate memory */
    freemem(&mpm, &mpw, &gv);

    if (gv.SEISMO)
        free_imatrix(recpos, 1, 3, 1, gv.NTRG);

    MPI_Barrier(MPI_COMM_WORLD);

    if (gv.MPID == 0) {
        time_av_v_update = time_av_v_update / (double)infocounter;
        time_av_s_update = time_av_s_update / (double)infocounter;
        time_av_v_exchange = time_av_v_exchange / (double)infocounter;
        time_av_s_exchange = time_av_s_exchange / (double)infocounter;
        time_av_timestep = time_av_timestep / (double)infocounter;
        log_info("Approximate average times for\n");
        log_info("  velocity update: .. %.6lfs.\n", time_av_v_update);
        log_info("  stress update: .... %.6lfs.\n", time_av_s_update);
        log_info("  velocity exchange:  %.6lfs.\n", time_av_v_exchange);
        log_info("  stress exchange: .. %.6lfs.\n", time_av_s_exchange);
        log_info("  time step: ........ %.6lfs.\n", time_av_timestep);
        cpu_time = clock() - cpu_time1;
        log_info("CPU time of program per PE: %.3lfs.\n", (double)cpu_time / (double)CLOCKS_PER_SEC);
        time8 = MPI_Wtime();
        log_info("Total real time of program: %.3lfs.\n", time8 - time1);
    }

    /* finalize logging */
    log_finalize();

    if (log_fp)
        fclose(log_fp);

    /* finalize MPI */
    MPI_Finalize();

    return EXIT_SUCCESS;
}
