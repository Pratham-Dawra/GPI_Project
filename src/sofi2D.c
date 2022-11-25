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
 *  Parallel 2-D Viscoelastic Finite Difference Seismic Modeling
 *  using the Standard Staggered Grid (SSG)
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
#ifdef EBUG
#include "debug_buffers.h"
#endif

#include <unistd.h>

int main ( int argc, char **argv )
{
    /* variables in main */
    int ns, nseismograms = 0, nt, nd, fdo3;
    int lsnap, nsnap = 0, lsamp = 0, buffsize;
    int ntr = 0, ntr_loc = 0, ntr_glob = 0, nsrc = 0, nsrc_loc = 0;
    int ishot, nshots; /* Added ishot and nshots for multiple shots */
    /*Limits for local grids defined in subgrid_bounds.c */
    char sigf[STRING_SIZE*2], file_ext[5];
    int * gx=NULL, * gy=NULL;
    clock_t cpu_time1 = 0, cpu_time = 0;
    FILE *log_fp = NULL;
    
    float memdyn, memmodel, memseismograms, membuffer, memtotal, memcpml=0.0;
    float fac1, fac2, memadd,memadd_L;
    char *buff_addr, ext[10];
    double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0;
    double time5 = 0.0, time6 = 0.0, time7 = 0.0, time8 = 0.0;
    double time_av_v_update = 0.0, time_av_s_update = 0.0,
    time_av_v_exchange = 0.0, time_av_s_exchange = 0.0, time_av_timestep = 0.0;
    
    float **psxx = NULL, **psxy = NULL, **psyy = NULL;
    float **pvx = NULL, **pvy = NULL, ***pr = NULL;
    float ***pp = NULL, ***pq = NULL;
    float **pu = NULL, **puipjp = NULL, **ptaus = NULL, **ptaup = NULL,
    *etaip = NULL, *etajm = NULL, *peta = NULL, **ptausipjp = NULL,
    **fipjp = NULL, ***dip = NULL, *bip = NULL, *bjm = NULL, *cip = NULL,
    *cjm = NULL, ***d = NULL, ***e = NULL, **f = NULL, **g = NULL;
    float **prho = NULL, **prip = NULL, **prjp = NULL, **ppi = NULL;
    
    float **pc11=NULL, **pc33=NULL, **pc13=NULL, **pc55=NULL, **pc15=NULL, **pc35=NULL;
    float **pc55ipjp=NULL, **pc15ipjp=NULL, **pc35ipjp=NULL;
    float **pc11u=NULL, **pc33u=NULL, **pc13u=NULL, **pc55u=NULL, **pc15u=NULL, **pc35u=NULL;
    float **pc55ipjpu=NULL, **pc15ipjpu=NULL, **pc35ipjpu=NULL;
    float ***pc11d=NULL, ***pc33d=NULL, ***pc13d=NULL, ***pc55d=NULL, ***pc15d=NULL, ***pc35d=NULL;
    float ***pc55ipjpd=NULL, ***pc15ipjpd=NULL, ***pc35ipjpd=NULL;

    float **ptau11=NULL, **ptau33=NULL, **ptau13=NULL, **ptau55=NULL, **ptau15=NULL, **ptau35=NULL;
    float **ptau55ipjp=NULL, **ptau15ipjp=NULL, **ptau35ipjp=NULL;
    
    float **pvxx=NULL, **pvyy=NULL, **pvyx=NULL, **pvxy=NULL;

    /* Save old spatial derivations of velocity for Adam Bashforth */
    float ** vxx_1=NULL,** vxx_2=NULL,** vxx_3=NULL,** vxx_4=NULL;
    float ** vyy_1=NULL,** vyy_2=NULL,** vyy_3=NULL,** vyy_4=NULL;
    float ** vxy_1=NULL,** vxy_2=NULL,** vxy_3=NULL,** vxy_4=NULL;
    float ** vyx_1=NULL,** vyx_2=NULL,** vyx_3=NULL,** vyx_4=NULL;
    
    /* Save old derivation of the stress for Adam Bashforth */
    float ** svx_1=NULL,** svx_2=NULL,** svx_3=NULL,** svx_4=NULL;
    float ** svy_1=NULL,** svy_2=NULL,** svy_3=NULL,** svy_4=NULL;

    /* Save old memory variables */
    float *** pr_2=NULL,*** pr_3=NULL,*** pr_4=NULL;
    float *** pq_2=NULL,*** pq_3=NULL,*** pq_4=NULL;
    float *** pp_2=NULL,*** pp_3=NULL,*** pp_4=NULL;
    
    /* We need some pointers for the time shift for Adam Bashforth*/
    float ** shift_s1=NULL,** shift_s2=NULL;
    float ** shift_v1=NULL,** shift_v2=NULL,** shift_v3=NULL,** shift_v4=NULL;
    float *** shift_r1=NULL,*** shift_r2=NULL,*** shift_r3=NULL;
    
    float **sectionvx = NULL, **sectionvy = NULL, **sectionp = NULL,
    **sectioncurl = NULL, **sectiondiv = NULL;
    float **absorb_coeff = NULL;
    float **srcpos = NULL, **srcpos_loc = NULL, **signals = NULL, *hc = NULL, **srcpos_current = NULL;
    int **recpos = NULL, **recpos_loc = NULL;
    
    float **bufferlef_to_rig = NULL, **bufferrig_to_lef = NULL,
    **buffertop_to_bot = NULL, **bufferbot_to_top = NULL;
    
    float ** seismo_fulldata=NULL;
    int * recswitch = NULL;
    
    /* PML variables */
    float * d_x=NULL, * K_x=NULL, * alpha_prime_x=NULL, * a_x=NULL, * b_x=NULL, * d_x_half=NULL,
    * K_x_half=NULL, * alpha_prime_x_half=NULL, * a_x_half=NULL, * b_x_half=NULL,
    * d_y=NULL, * K_y=NULL, * alpha_prime_y=NULL, * a_y=NULL, * b_y=NULL, * d_y_half=NULL,
    * K_y_half=NULL, * alpha_prime_y_half=NULL, * a_y_half=NULL, * b_y_half=NULL;
    float ** psi_sxx_x=NULL, ** psi_syy_y=NULL, ** psi_sxy_y=NULL, ** psi_sxy_x=NULL,
    ** psi_vxx=NULL, ** psi_vyy=NULL, ** psi_vxy=NULL, ** psi_vyx=NULL, ** psi_vxxs=NULL;

    MPI_Request *req_send, *req_rec;
    
    GlobVar gv = {.MPID=-1, .OUTNTIMESTEPINFO=100, .NDT=1, .IDX=1, .IDY=1};  
    
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
	log_fatal("Unexpected number of commandline arguments; single argument required: name of json parameter file.\n");
      }

      const char *fileinp = argv[1];

      /* check if parameter file can be opened */
      if (access(fileinp, R_OK) != 0) {
	log_fatal("Cannot open/read json parameter file %s.\n", fileinp);
      }
      
      /* check suffix of parameter file */
      if (!STRSTRCOMP(fileinp,".json")) {
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
    case 0:   /* logging to stdout/stderr for all ranks */
      log_fp = NULL;
      log_set_output(NULL);
      log_infoc(0, "Log messages sent to stdout/stderr on all MPI ranks.\n");
      break;
    case 1:  /* logging to file for all ranks */
      log_infoc(0, "Now redirecting log messages to log file on all MPI ranks.\n");
      if ((log_fp=fopen(gv.LOG_FILE, "w")) == NULL) {
	log_fatal("Opening log file %s for writing failed.\n", gv.LOG_FILE);
      }
      log_set_output(log_fp);
      log_info("This is the log file %s generated by PE %d.\n", gv.LOG_FILE, gv.MPID);
      break;
    case 2:  /* logging to stdout/stderr on rank 0, logging to file for all other ranks */
      if (0==gv.MPID) {
	log_fp = NULL;
	log_set_output(NULL);
	log_info("Now redirecting log messages to log file on all MPI ranks except rank 0.\n");
      } else {
	if ((log_fp=fopen(gv.LOG_FILE, "w")) == NULL) {
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

    gv.NT = iround(gv.TIME/gv.DT);      /* number of time steps */
    ns = iround(gv.NT/gv.NDT);          /* number of samples per trace */
    lsnap = iround(gv.TSNAP1/gv.DT);    /* first snapshot at this time step */
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
    
    if (gv.SEISMO) {
        recpos = receiver(&ntr, &gv);
        recswitch = ivector(1, ntr);
        recpos_loc = splitrec(recpos, &ntr_loc, ntr, recswitch, &gv);
        ntr_glob = ntr;
        ntr = ntr_loc;
    }
    
    /* allocate buffer for seismogram output, merged seismogram section of all PEs */
    if ( gv.SEISMO ) seismo_fulldata=matrix ( 1,ntr_glob,1,ns );
    
    /* number of seismogram sections which have to be stored in core memory*/
    /* allocation of memory for seismogramm merge */
    switch ( gv.SEISMO ) {
        case 1: /* particle velocities only */
            nseismograms = 3;
            break;
        case 2: /* pressure only */
            nseismograms = 1;
            break;
        case 3: /* curl and div only */
            nseismograms = 2;
            break;
        case 4: /* everything */
            nseismograms = 6;
            break;
        default:
            nseismograms = 1;
            break;
    }
    
    /*allocate memory for dynamic, static and buffer arrays */
    nd = gv.FDORDER/2;
    fdo3 = 2 * nd;
    
    fac1 = ( gv.NX + fdo3 ) * ( gv.NY + fdo3 );
    fac2 = sizeof ( float ) * pow ( 2.0, -20.0 );
    memadd=0.0; memadd_L=0.0;
    if(gv.FDORDER_TIME==4){memadd=24.0; memadd_L=9;}
    if ( gv.L ) {
        memdyn = ( memadd+ 5.0 + (3.0+memadd_L) * ( float ) gv.L ) * fac1 * fac2;
        memmodel = ( 12.0 + 3.0 * ( float ) gv.L ) * fac1 * fac2 + gv.NX * gv.NY * fac2;
    } else {
        memdyn = (memadd+5.0) * fac1 * fac2;
        memmodel = 6.0 * fac1 * fac2 + gv.NX * gv.NY * fac2;
    }
    
    memseismograms = nseismograms * ntr * ns * fac2;
    membuffer = 2.0 * fdo3 * ( gv.NY + gv.NX ) * fac2;
    buffsize = 2.0 * 2.0 * fdo3 * ( gv.NX + gv.NY ) * sizeof ( MPI_FLOAT );
    if ( gv.ABS_TYPE==1 ) memcpml=2.0*gv.FW*4.0* ( gv.NY+gv.NX ) *fac2+20.0*2.0*gv.FW*fac2;
    memtotal = memdyn + memmodel + memseismograms + membuffer +memcpml
    + ( buffsize * pow ( 2.0, -20.0 ) );
    
   if (gv.MPID == 0) {
      log_info("Size of local grids: NX=%d, NY=%d\n", gv.NX, gv.NY);
      log_info("Each process is now trying to allocate memory for:\n");
      log_info("  Dynamic variables: ............. %6.2f MB\n", memdyn);
      log_info("  Static variables: .............. %6.2f MB\n", memmodel);
      log_info("  Seismograms: ................... %6.2f MB\n", memseismograms);
      log_info("  Buffer arrays for grid exchange: %6.2f MB\n", membuffer);
      log_info("  Network buffer for MPI_Bsend: .. %6.2f MB\n", buffsize*pow(2.0,-20.0));
      if (gv.ABS_TYPE == 1) log_info("  CPML variables: ................ %6.2f MB\n", memcpml);
      log_info("------------------------------------------------\n");
      log_info("Total memory required: ........... %6.2f MB.\n", memtotal);
    }

    /* allocate buffer for buffering messages */
    buff_addr = malloc ( buffsize );
    if ( !buff_addr ) log_fatal("Allocation failure for MPI_Bsend buffer!\n");
    MPI_Buffer_attach ( buff_addr, buffsize );
    
    /* allocation for request and status arrays */
    req_send = ( MPI_Request * ) malloc ( REQUEST_COUNT * sizeof ( MPI_Request ) );
    req_rec = ( MPI_Request * ) malloc ( REQUEST_COUNT * sizeof ( MPI_Request ) );
    /*	send_statuses = (MPI_Status *) malloc(REQUEST_COUNT * sizeof(MPI_Status));
     rec_statuses = (MPI_Status *) malloc(REQUEST_COUNT * sizeof(MPI_Status)); */
    
    /* ------------ memory allocation for arrays ------------- */
    /* subgrid arrays*/
    gy = ivector ( 1,4 );
    gx = ivector ( 1,4 );
    
    /* dynamic (wavefield) arrays (elastic + viscoelastic) */    
    psxx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    psxy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    psyy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    pvx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    pvy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    
    if (gv.FDORDER_TIME==4) {
        vxx_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxx_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxx_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxx_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        vyy_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyy_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyy_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyy_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        vxy_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxy_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxy_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vxy_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        vyx_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyx_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyx_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        vyx_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        svx_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svx_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svx_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svx_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        svy_1= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svy_2= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svy_3= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        svy_4= matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    }
    
    if ( gv.ABS_TYPE==1 ) {
        /* PML */
        d_x = vector ( 1,2*gv.FW );
        K_x = vector ( 1,2*gv.FW );
        alpha_prime_x = vector ( 1,2*gv.FW );
        a_x = vector ( 1,2*gv.FW );
        b_x = vector ( 1,2*gv.FW );
        
        d_x_half = vector ( 1,2*gv.FW );
        K_x_half = vector ( 1,2*gv.FW );
        alpha_prime_x_half = vector ( 1,2*gv.FW );
        a_x_half = vector ( 1,2*gv.FW );
        b_x_half = vector ( 1,2*gv.FW );
        
        d_y = vector ( 1,2*gv.FW );
        K_y = vector ( 1,2*gv.FW );
        alpha_prime_y = vector ( 1,2*gv.FW );
        a_y = vector ( 1,2*gv.FW );
        b_y = vector ( 1,2*gv.FW );
        
        d_y_half = vector ( 1,2*gv.FW );
        K_y_half = vector ( 1,2*gv.FW );
        alpha_prime_y_half = vector ( 1,2*gv.FW );
        a_y_half = vector ( 1,2*gv.FW );
        b_y_half = vector ( 1,2*gv.FW );
        
        psi_sxx_x =  matrix ( 1,gv.NY,1,2*gv.FW );
        psi_syy_y =  matrix ( 1,2*gv.FW,1,gv.NX );
        psi_sxy_y =  matrix ( 1,2*gv.FW,1,gv.NX );
        psi_sxy_x =  matrix ( 1,gv.NY,1,2*gv.FW );
        
        psi_vxx   =  matrix ( 1,gv.NY,1,2*gv.FW );
        psi_vyy   =  matrix ( 1,2*gv.FW,1,gv.NX );
        psi_vxy   =  matrix ( 1,2*gv.FW,1,gv.NX );
        psi_vyx   =  matrix ( 1,gv.NY,1,2*gv.FW );
        
        psi_vxxs  =  matrix ( 1,gv.NY,1,2*gv.FW ); /* For surface_elastic(visc).c*/
    }
    
    /* dynamic (wavefield) arrays (viscoelastic) */
    if ( gv.L > 0 ) {
        pr = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pp = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pq = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
    }
    
    if(gv.L>0 && gv.FDORDER_TIME==4){
        pr_2 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pr_3 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pr_4 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        
        pp_2 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pp_3 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pp_4 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        
        pq_2 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pq_3 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pq_4 = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
    }
    
    /* static (model) arrays (isotropic elastic + viscoelastic) */
    if (gv.WEQ>2){
        prho = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        prip = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        prjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ppi = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pu = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        puipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        absorb_coeff = matrix ( 1, gv.NY, 1, gv.NX );
    }
    
    /* static (model) arrays (viscoelastic) */
    if ( gv.WEQ==4 ) { /*viscoelastic  isotropic wave equation */
        dip = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        e = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        ptaus = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptausipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptaup = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        fipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        f = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        g = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        peta = vector ( 1, gv.L );
        etaip = vector ( 1, gv.L );
        etajm = vector ( 1, gv.L );
        bip = vector ( 1, gv.L );
        bjm = vector ( 1, gv.L );
        cip = vector ( 1, gv.L );
        cjm = vector ( 1, gv.L );
    }
    
    if ( gv.WEQ==5 ) {/*elastic VTI wave equation */
        pc11 = ppi;
        pc33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc55 = pu;
        pc55ipjp = puipjp;
    }   

    if ( gv.WEQ==6 ) { /*viscoelastic VTI wave equation */
        pc11 = ppi;
        pc33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc55 = pu;
        pc55ipjp = puipjp;
 
        ptau11 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau15 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau55 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau55ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        pc55ipjpd = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc13d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc33d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc11d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        
        pc55ipjpu = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc11u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc33u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        peta = vector ( 1, gv.L );
        bip = vector ( 1, gv.L );
        cip = vector ( 1, gv.L );

    }
    
    if ( gv.WEQ==7 ) {/*elastic TTI wave equation */
        pc11 = ppi;
        pc33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc55 = pu;
        pc55ipjp = puipjp;
        pc15 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc15ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

        pvxx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvyy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvyx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvxy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    }
    
    if ( gv.WEQ==8 ) { /*viscoelastic TTI wave equation */
        pc11 = ppi;
        pc33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc55 = pu;
        pc15 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        pc55ipjp = puipjp;
        pc15ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

 
        ptau11 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau33 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau13 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau55 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau15 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau35 = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau55ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau15ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        ptau35ipjp = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
  

        pc11d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc33d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc13d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc55d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc15d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc35d = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc55ipjpd = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc15ipjpd = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        pc35ipjpd = f3tensor ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
 
        pc11u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc33u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc13u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc55u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc15u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35u = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        pc55ipjpu = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc15ipjpu = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pc35ipjpu = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
 
        peta = vector ( 1, gv.L );
        bip = vector ( 1, gv.L );
        cip = vector ( 1, gv.L );
        
        pvxx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvyy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvyx = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        pvxy = matrix ( -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

    }
    
    /* memory allocation for buffer arrays in which the wavefield
     information to be exchanged between neighboring PEs is stored */
        
    bufferlef_to_rig = matrix ( 1, gv.NY, 1, fdo3 );
    bufferrig_to_lef = matrix ( 1, gv.NY, 1, fdo3 );
    buffertop_to_bot = matrix ( 1, gv.NX, 1, fdo3 );
    bufferbot_to_top = matrix ( 1, gv.NX, 1, fdo3 );
      
    if ( ntr > 0 ) {
        switch ( gv.SEISMO ) {
            case 1: /* particle velocities only */
                sectionvx = matrix ( 1, ntr, 1, ns );
                sectionvy = matrix ( 1, ntr, 1, ns );
                break;
            case 2: /* pressure only */
                sectionp = matrix ( 1, ntr, 1, ns );
                break;
            case 3: /* curl and div only */
                sectioncurl = matrix ( 1, ntr, 1, ns );
                sectiondiv = matrix ( 1, ntr, 1, ns );
                break;
            case 4: /* everything */
                sectionvx = matrix ( 1, ntr, 1, ns );
                sectionvy = matrix ( 1, ntr, 1, ns );
                sectioncurl = matrix ( 1, ntr, 1, ns );
                sectiondiv = matrix ( 1, ntr, 1, ns );
                sectionp = matrix ( 1, ntr, 1, ns );
                break;
        }
    }
    
    /* memory for source position definition for saving the current positions */
    srcpos_current = matrix ( 1, NSPAR, 1, 1 );
    
    log_debug("Memory allocation for PE %d was successful.\n", gv.MPID);
    
    /* Holberg coefficients for FD operators*/
    hc = holbergcoeff(&gv);
    /* Reading source positions from SOURCE_FILE */
    srcpos = sources(&nsrc, &gv);

    MPI_Barrier(MPI_COMM_WORLD);

    /* create model grids */
    if ( gv.READMOD ){
        switch ( gv.WEQ ) {
            case 3: /* elastic */
                readmod_elastic(prho, ppi, pu, &gv);
                break;
            case 4: /* viscoelastic */
                readmod_visco(prho, ppi, pu, ptaus, ptaup, peta, &gv);
                break;
            case 5: /* elastic VTI */
              readmod_elastic_vti(prho, pc11, pc33, pc13, pc55, &gv);
#ifdef EBUG
	      debug_check_matrix(prho, 0, gv.NX, gv.NY, 5, 0, "prho");
	      debug_check_matrix(pc11, 0, gv.NX, gv.NY, 5, 0, "pc11");
	      debug_check_matrix(pc33, 0, gv.NX, gv.NY, 5, 0, "pc33");
	      debug_check_matrix(pc13, 0, gv.NX, gv.NY, 5, 0, "pc13");
	      debug_check_matrix(pc55, 0, gv.NX, gv.NY, 5, 0, "pc55");
#endif
                break;
            case 6: /* viscoelastic VTI */
              readmod_visco_vti(prho, pc11, pc33, pc13, pc55, ptau11, ptau33, ptau13, ptau55, peta, &gv);
#ifdef EBUG
	      debug_check_matrix(prho, 0, gv.NX, gv.NY, 6, 0, "prho");
	      debug_check_matrix(pc11, 0, gv.NX, gv.NY, 6, 0, "pc11");
	      debug_check_matrix(pc33, 0, gv.NX, gv.NY, 6, 0, "pc33");
	      debug_check_matrix(pc13, 0, gv.NX, gv.NY, 6, 0, "pc13");
	      debug_check_matrix(pc55, 0, gv.NX, gv.NY, 6, 0, "pc55");
	      debug_check_matrix(ptau11, 0, gv.NX, gv.NY, 6, 0, "ptau11");
	      debug_check_matrix(ptau33, 0, gv.NX, gv.NY, 6, 0, "ptau33");
	      debug_check_matrix(ptau13, 0, gv.NX, gv.NY, 6, 0, "ptau13");
	      debug_check_matrix(ptau55, 0, gv.NX, gv.NY, 6, 0, "ptau55");
	      debug_check_vector(peta, 0, gv.L, 6, 0, "peta");
#endif
                break;
            case 7: /* elastic TTI */
              readmod_elastic_tti(prho, pc11, pc33, pc13, pc55, pc15, pc35, &gv);
#ifdef EBUG
	      debug_check_matrix(prho, 0, gv.NX, gv.NY, 7, 0, "prho");
	      debug_check_matrix(pc11, 0, gv.NX, gv.NY, 7, 0, "pc11");
	      debug_check_matrix(pc33, 0, gv.NX, gv.NY, 7, 0, "pc33");
	      debug_check_matrix(pc13, 0, gv.NX, gv.NY, 7, 0, "pc13");
	      debug_check_matrix(pc55, 0, gv.NX, gv.NY, 7, 0, "pc55");
	      debug_check_matrix(pc15, 0, gv.NX, gv.NY, 7, 0, "pc15");
	      debug_check_matrix(pc35, 0, gv.NX, gv.NY, 7, 0, "pc35");
#endif
                break;
            case 8: /* viscoelastic TTI */
              readmod_visco_tti ( prho, pc11, pc33, pc13, pc55, pc15, pc35,
                                 ptau11, ptau33, ptau13, ptau55, ptau15, ptau35, peta, &gv );
#ifdef EBUG
	      debug_check_matrix(prho, 0, gv.NX, gv.NY, 8, 0, "prho");
	      debug_check_matrix(pc11, 0, gv.NX, gv.NY, 8, 0, "pc11");
	      debug_check_matrix(pc33, 0, gv.NX, gv.NY, 8, 0, "pc33");
	      debug_check_matrix(pc13, 0, gv.NX, gv.NY, 8, 0, "pc13");
	      debug_check_matrix(pc55, 0, gv.NX, gv.NY, 8, 0, "pc55");
	      debug_check_matrix(pc15, 0, gv.NX, gv.NY, 8, 0, "pc15");
	      debug_check_matrix(pc35, 0, gv.NX, gv.NY, 8, 0, "pc35");
	      debug_check_matrix(ptau11, 0, gv.NX, gv.NY, 8, 0, "ptau11");
	      debug_check_matrix(ptau33, 0, gv.NX, gv.NY, 8, 0, "ptau33");
	      debug_check_matrix(ptau13, 0, gv.NX, gv.NY, 8, 0, "ptau13");
	      debug_check_matrix(ptau55, 0, gv.NX, gv.NY, 8, 0, "ptau55");
	      debug_check_matrix(ptau15, 0, gv.NX, gv.NY, 8, 0, "ptau15");
	      debug_check_matrix(ptau35, 0, gv.NX, gv.NY, 8, 0, "ptau35");
	      debug_check_vector(peta, 0, gv.L, 8, 0, "peta");
#endif
	      break;
        }
    }else {
        switch ( gv.WEQ ) {
            case 3: /* elastic */
                model_elastic ( prho, ppi, pu, &gv );
                break;
            case 4: /* viscoelastic */
                model_visco ( prho, ppi, pu, ptaus, ptaup, peta, &gv );
                break;
            case 5 : /* elastic VTI */
               model_elastic_VTI ( prho, pc11, pc33, pc13, pc55, &gv );
                break;

            case 6 : /* viscoelastic VTI */
               model_visco_vti ( prho, pc11, pc33, pc13, pc55, ptau11, ptau33, ptau13, ptau55, peta, &gv );
                break;
                
            case 7 : /* elastic TTI */
               model_elastic_TTI ( prho, pc11, pc33, pc13, pc55, pc15, pc35, &gv );
                break;

            case 8 : /* viscoelastic TTI */
               model_visco_tti ( prho, pc11, pc33, pc13, pc55, pc15, pc35, ptau11,
                                ptau33, ptau13, ptau55, ptau15, ptau35, peta, &gv );
                break;
        }
    }

    MPI_Barrier ( MPI_COMM_WORLD );

    /* check if the FD run will be stable and free of numerical dispersion */
    switch ( gv.WEQ) {
            case 3: /* elastic */
	        checkfd (prho, ppi, pu, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
                break;
            case 4: /* viscoelastic */
                checkfd (prho, ppi, pu, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
                break;
            case 5 : /* elastic VTI */
                checkfd (prho, pc11, pc55, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
                break;
             case 6 : /* viscoelastic VTI */
                checkfd (prho, pc11, pc55, ptau55, ptau11, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
                break;
            case 7 : /* elastic TTI */
                checkfd (prho, pc11, pc55, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
		break;
            case 8 : /* viscoelastic TTI */
                checkfd (prho, pc11, pc55, ptau55, ptau11, peta, hc, srcpos, nsrc, recpos, ntr_glob, &gv );
		break;
        }
    
    /* calculate damping coefficients for CPMLs*/
    if ( gv.ABS_TYPE==1 ) {
        PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half,
		d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half, &gv);
    }

    MPI_Barrier(MPI_COMM_WORLD );
    
    /* calculate 2-D array for exponential damping of reflections at the edges of the numerical mesh */
    if ( gv.ABS_TYPE==2 ) {
      absorb(absorb_coeff, &gv);
    }
    
    /* For the calculation of the material parameters beteween gridpoints
     they have to be averaged. For this, values lying at 0 and NX+1,
     for example, are required on the local grid. These are now copied from the
     neighbouring grids */

    switch ( gv.WEQ ) {
        case 3: /* elastic */
            matcopy_elastic (prho, ppi, pu, &gv );
            av_mue ( pu, puipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
            break;
        case 4: /* viscoelastic */
            matcopy ( prho, ppi, pu, ptaus, ptaup, &gv );
            av_mue ( pu, puipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
            av_tau ( ptaus, ptausipjp, &gv );
            break;
        case 5 : /* elastic VTI */
            matcopy_elastic ( prho, pc11, pc55, &gv );
            av_mue ( pc55, pc55ipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
#ifdef EBUG
	    debug_check_matrix(pc55ipjp, 0, gv.NX, gv.NY, 55, 0, "pc55ipjp");
	    debug_check_matrix(prip, 0, gv.NX, gv.NY, 55, 0, "prip");
	    debug_check_matrix(prjp, 0, gv.NX, gv.NY, 55, 0, "prjp");
#endif
            break;
        case 6 : /* viscoelastic VTI */
            matcopy_elastic ( prho, ptau55, pc55, &gv );
            av_mue ( pc55, pc55ipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
            av_tau (ptau55, ptau55ipjp, &gv );
#ifdef EBUG
	    debug_check_matrix(pc55ipjp, 0, gv.NX, gv.NY, 66, 0, "pc55ipjp");
	    debug_check_matrix(prip, 0, gv.NX, gv.NY, 66, 0, "prip");
	    debug_check_matrix(prjp, 0, gv.NX, gv.NY, 66, 0, "prjp");
	    debug_check_matrix(ptau55ipjp, 0, gv.NX, gv.NY, 66, 0, "ptau55ipjp");
#endif
            break;
        case 7 : /* elastic TTI */
            matcopy_elastic ( prho, pc11, pc55, &gv );
            matcopy_elastic ( prho, pc15, pc35, &gv );
            av_mue ( pc55, pc55ipjp, &gv );
            av_mue ( pc15, pc15ipjp, &gv );
            av_mue ( pc35, pc35ipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
#ifdef EBUG
	    debug_check_matrix(pc55ipjp, 0, gv.NX, gv.NY, 77, 0, "pc55ipjp");
	    debug_check_matrix(prip, 0, gv.NX, gv.NY, 77, 0, "prip");
	    debug_check_matrix(prjp, 0, gv.NX, gv.NY, 77, 0, "prjp");
	    debug_check_matrix(pc15ipjp, 0, gv.NX, gv.NY, 77, 0, "pc15ipjp");
	    debug_check_matrix(pc35ipjp, 0, gv.NX, gv.NY, 77, 0, "pc35ipjp");
#endif
            break;
        case 8 : /* viscoelastic TTI */
            matcopy_elastic ( prho, ptau55, pc55, &gv );
            matcopy_elastic ( prho, ptau15, pc15, &gv );
            matcopy_elastic ( prho, ptau35, pc35, &gv );
            av_mue ( pc55, pc55ipjp, &gv );
            av_mue ( pc15, pc15ipjp, &gv );
            av_mue ( pc35, pc35ipjp, &gv );
            av_rho ( prho, prip, prjp, &gv );
            av_tau (ptau55, ptau55ipjp, &gv );
            av_tau (ptau15, ptau15ipjp, &gv );
            av_tau (ptau35, ptau35ipjp, &gv );
#ifdef EBUG
	    debug_check_matrix(pc55ipjp, 0, gv.NX, gv.NY, 88, 0, "pc55ipjp");
	    debug_check_matrix(prip, 0, gv.NX, gv.NY, 88, 0, "prip");
	    debug_check_matrix(prjp, 0, gv.NX, gv.NY, 88, 0, "prjp");
	    debug_check_matrix(pc15ipjp, 0, gv.NX, gv.NY, 88, 0, "pc15ipjp");
	    debug_check_matrix(pc35ipjp, 0, gv.NX, gv.NY, 88, 0, "pc35ipjp");
	    debug_check_matrix(ptau55ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau55ipjp");
	    debug_check_matrix(ptau15ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau15ipjp");
	    debug_check_matrix(ptau35ipjp, 0, gv.NX, gv.NY, 88, 0, "ptau35ipjp");
#endif
	    break;
    }

    
    /* Preparing memory variables for update_s (viscoelastic only) */

    if (gv.FDORDER_TIME==2) {
        switch ( gv.WEQ ) {
            case 4:  /* viscoelastic */
                prepare_update_s ( etajm, etaip, peta, fipjp, pu, puipjp, ppi, ptaus,
                                  ptaup, ptausipjp, f, g, bip, bjm, cip, cjm, dip, d, e, &gv );
                break;
            case 5:
                dt_mult(gv.NX,gv.NY,gv.DT,pc11);
                dt_mult(gv.NX,gv.NY,gv.DT,pc33);
                dt_mult(gv.NX,gv.NY,gv.DT,pc13);
                dt_mult(gv.NX,gv.NY,gv.DT,pc55ipjp);
                break;
            case 6:  /* viscoelastic VTI*/
                
                prepare_update_s_vti(peta, pc11, pc13, pc33,  pc55ipjp,
                                     ptau11, ptau13, ptau33,  ptau55ipjp,
                                     pc55ipjpu, pc13u, pc11u,  pc33u,
                                     pc55ipjpd,  pc13d, pc11d,  pc33d,
                                     bip,  cip, &gv);
                
                break;
                
            case 7:
                dt_mult(gv.NX,gv.NY,gv.DT,pc11);
                dt_mult(gv.NX,gv.NY,gv.DT,pc33);
                dt_mult(gv.NX,gv.NY,gv.DT,pc13);
                dt_mult(gv.NX,gv.NY,gv.DT,pc15);
                dt_mult(gv.NX,gv.NY,gv.DT,pc35);
                dt_mult(gv.NX,gv.NY,gv.DT,pc55ipjp);
                dt_mult(gv.NX,gv.NY,gv.DT,pc35ipjp);
                dt_mult(gv.NX,gv.NY,gv.DT,pc15ipjp);
                break;
                
            case 8:  /* viscoelastic TTI*/
                
                prepare_update_s_tti(peta, pc11, pc33, pc13, pc55, pc15, pc35, pc55ipjp, pc15ipjp, pc35ipjp,
                                     ptau11, ptau33, ptau13, ptau55, ptau15, ptau35, ptau55ipjp, ptau15ipjp, ptau35ipjp,
                                     pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                     pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                     bip,  cip, &gv);
                
                break;
            }
    }

    if ( (gv.WEQ==4) && gv.FDORDER_TIME==4) {
        prepare_update_s_4 ( etajm, etaip, peta, fipjp, pu, puipjp, ppi, ptaus,
                          ptaup, ptausipjp, f, g, bip, bjm, cip, cjm, dip, d, e, &gv );
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* comunication initialisation for persistent communication */
    /*comm_ini(bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
     req_send, req_rec);*/
    /* currently MPI_Sendrecv_replace is used! */
    
    time2 = MPI_Wtime();
    log_infoc(0, "Starting time stepping around real time %4.2fs.\n", time2-time1);
    
    /*----------------------  loop over multiple shots  ------------------*/
    
    if ( gv.RUN_MULTIPLE_SHOTS ) {
        nshots = nsrc;
    } else {
        nshots = 1;
    }

    for ( ishot = 1; ishot <= nshots; ishot++ ) {
        
        for ( nt = 1; nt <= 12; nt++ ) {
            srcpos_current[nt][1] = srcpos[nt][ishot];
        }
        
        if ( gv.RUN_MULTIPLE_SHOTS ) {
	  log_info("Starting simulation for shot %d of %d.\n", ishot, nshots);
	  //log_info("number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\n");
	  //log_info("   %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f  \t %6.2f\n", ishot, srcpos_current[1][1], srcpos_current[2][1],
	  //         srcpos_current[4][1], srcpos_current[5][1], srcpos_current[6][1], srcpos_current[7][1]);
            
	  /* find this single source positions on subdomains  */
	  if ( nsrc_loc > 0 )
	    free_matrix ( srcpos_loc, 1, 12, 1, 1 );
	  srcpos_loc = splitsrc ( srcpos_current, &nsrc_loc, 1, &gv );
        }
        else {
	  srcpos_loc = splitsrc ( srcpos, &nsrc_loc, nsrc, &gv ); /* Distribute source positions on subdomains */
        }
	
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* calculate wavelet for each source point */
        signals = wavelet ( srcpos_loc, nsrc_loc, &gv );

        /* write source wavelet to file */
        if (0==gv.MPID && 1==gv.SIGOUT) {
            switch (gv.SIGOUT_FORMAT){
                case 1: sprintf(file_ext,"su");  break;
                case 2: sprintf(file_ext,"txt"); break;
                case 3: sprintf(file_ext,"bin"); break;
            }
            sprintf(sigf,"%s.shot%d.%s",gv.SIGOUT_FILE,ishot,file_ext);
            log_info("Writing source wavelet to file %s.\n",sigf);
            outseis_glob(fopen(sigf,"w"),signals,recpos,recpos_loc,1,srcpos_loc,nsrc,gv.NT,gv.SIGOUT_FORMAT,ishot,0, &gv);
        }

        /* initialize wavefield with zero */
        
        if ( gv.ABS_TYPE == 1 ) {
            if ( gv.L )
                zero_PML_visc ( -nd+1,gv.NY+nd,-nd+1,gv.NX+nd,pvx,pvy,psxx,psyy,psxy,psi_sxx_x,psi_sxy_x,
                               psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pr,pp,pq, &gv );
            else
                zero_PML_elastic ( -nd+1,gv.NY+nd,-nd+1,gv.NX+nd,pvx,pvy,psxx,psyy,psxy,psi_sxx_x,psi_sxy_x,
                                  psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs, &gv );
        }
        
        if(gv.FDORDER_TIME==4){
            if(gv.L) {
                zero_visco_4(-nd+1,gv.NY+nd,-nd+1,gv.NX+nd,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,
				vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,
				svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4, &gv );
            } else {
                zero_elastic_4(-nd+1,gv.NY+nd,-nd+1,gv.NX+nd,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,
				vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4 );
            }
        }
        
        if ( gv.ABS_TYPE != 1 ) {
            if ( gv.L )
                zero_visc ( -nd+1, gv.NX+nd, -nd+1,gv.NY+nd, pvx, pvy, psxx, psyy, psxy, pr, pp, pq, &gv );
            else
                zero_elastic ( -nd+1, gv.NX + nd, -nd+1, gv.NY+nd, pvx, pvy, psxx, psyy, psxy );
        }
        
        /* Reseting lsmap to NDT for saving seismograms  */
        lsamp = gv.NDT;
        
        subgrid_bounds ( 1, gv.NX, 1, gv.NY, gx, gy, &gv );
        /*---------------------------------------------------------------*/
        /*----------------------  loop over timesteps  ------------------*/
        /*---------------------------------------------------------------*/
        
        for ( nt = 1; nt <= gv.NT; nt++ ) {
            
            if (isnan(pvy[gv.NY/2][gv.NX/2])) {
	      log_error("Time step: %d; pvy: %f.\n",nt,pvy[gv.NY/2][gv.NX/2]);
	      log_fatal("Simulation is unstable!\n");
            }

            if ( ( gv.MPID == 0 ) && ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
	      log_info("Computing time step %d of %d.\n", nt, gv.NT);
	      time3 = MPI_Wtime();
            }
            /*---------------------------------------------------------------*/
            /* update of particle velocities --------------------------------*/
            /*---------------------------------------------------------------*/
            if (gv.FDORDER_TIME==2) {
                update_v_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prho, prip, prjp,
                                   srcpos_loc, signals, nsrc_loc, hc, &gv );
#ifdef EBUG
		debug_check_matrix(pvx, nt, gv.NX, gv.NY, 121, 0, "pvx");
		debug_check_matrix(pvy, nt, gv.NX, gv.NY, 121, 0, "pvy");
#endif
                
                if ( gv.FW ) {
                    if ( gv.ABS_TYPE==1 ) {
                        update_v_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prip, prjp, hc,
                                      K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half,
                                      a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x, &gv );
                    }

                    if ( gv.ABS_TYPE == 2 ) {
                        update_v_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prip, prjp, absorb_coeff, hc, &gv );
                    }
#ifdef EBUG
		    debug_check_matrix(pvx, nt, gv.NX, gv.NY, 122, 0, "pvx");
		    debug_check_matrix(pvy, nt, gv.NX, gv.NY, 122, 0, "pvy");
#endif
                }
            }

            if (gv.FDORDER_TIME==4) {
                update_v_interior_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prho, prip, prjp,
                                        srcpos_loc, signals, nsrc_loc, hc,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4, &gv );
                if ( gv.FW ) {
                    if ( gv.ABS_TYPE==1 ) {
                        update_v_PML_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prip, prjp, hc,
                                           K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half,
                                           a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4, &gv );
                    }

                    if ( gv.ABS_TYPE == 2 ) {
                        update_v_abs_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, prip, prjp, absorb_coeff, hc,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4, &gv );
                    }
                }
                
                /* Shift spatial derivations of the stress one time step back */
                shift_s1=svx_4;svx_4=svx_3;svx_3=svx_2;svx_2=svx_1;svx_1=shift_s1;
                shift_s2=svy_4;svy_4=svy_3;svy_3=svy_2;svy_2=svy_1;svy_1=shift_s2;
            }
        
            if ( ( gv.MPID == 0 )
                && ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
                time4 = MPI_Wtime();
                time_av_v_update += ( time4 - time3 );
                log_debug("Starting particle velocity exchange between PEs...\n");
            }
            
            /*---------------------------------------------------------------*/
            /* ------- exchange of particle velocities between PEs --------------*/
            /*---------------------------------------------------------------*/

            exchange_v ( nd, pvx, pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec, &gv );
            
            if ((gv.WEQ==7) || (gv.WEQ==8)) { /* TTI */
                v_derivatives(pvx, pvy, pvxx, pvyy, pvyx, pvxy,hc, &gv);
                exchange_v ( nd, pvxx, pvyy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec, &gv );
                exchange_v ( nd, pvyx, pvxy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec, &gv );
            }
            
            if ( ( gv.MPID == 0 )
                && ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
                time5 = MPI_Wtime();
                time_av_v_exchange += ( time5 - time4 );
                log_debug("Finished particle velocity exchange between PEs (real time: %4.3fs).\n", time5 - time4);
            }
            /*---------------------------------------------------------------*/
            /* stress update ------------------------------------------------*/
            /*---------------------------------------------------------------*/
            if(gv.FDORDER_TIME==2){

                switch ( gv.WEQ) {
                     case 3: /* elastic */
                        update_s_elastic_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, ppi, pu, puipjp, hc, &gv );
                    
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_elastic_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, ppi, pu, puipjp, hc,
                                                  K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                        if ( gv.ABS_TYPE ==2 )
                            update_s_elastic_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy,
                                                  ppi, pu, puipjp, absorb_coeff, hc, &gv );
                    }
                    break;

                case 4:  /* viscoelastic */
                    update_s_visc_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pr, pp, pq, fipjp,
                                            f, g, bip, bjm, cip, cjm, d, e, dip, hc, &gv );
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_visc_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, hc, pr, pp, pq, fipjp,
                                               f, g, bip, bjm, cip, cjm, d, e, dip,K_x, a_x, b_x, K_x_half, a_x_half, b_x_half,
                                               K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                        if ( gv.ABS_TYPE ==2 )
                            update_s_visc_abs ( 1, gv.NX, 1, gv.NY, gx,gy, nt, pvx, pvy, psxx, psyy, psxy, pr,
                                               pp, pq, ppi, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
                                               absorb_coeff,hc, &gv );
                    }
                    break;
                   case 5: /* elastic VTI */
                        update_s_elastic_VTI_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy,
                                                       pc11, pc55ipjp, pc13, pc33, hc, &gv );
                    
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_elastic_VTI_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pc11, pc13, pc33, pc55ipjp, hc,
                                                  K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                        if ( gv.ABS_TYPE ==2 )
                            update_s_elastic_vti_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy,
                                                  pc11, pc55ipjp, pc13, pc33, absorb_coeff, hc, &gv );
                    }
#ifdef EBUG
		    debug_check_matrix(psxx, nt, gv.NX, gv.NY, 555, 0, "psxx");
		    debug_check_matrix(psyy, nt, gv.NX, gv.NY, 555, 0, "psyy");
		    debug_check_matrix(psxy, nt, gv.NX, gv.NY, 555, 0, "psxy");
#endif
                    break;

                    case 6: /* viscoelastic VTI */
                        update_s_visc_VTI_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pr, pp, pq,
                                                            pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
                                                            bip,  cip, hc, &gv);
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_visc_VTI_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pr, pp, pq,
                                                            pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
                                                            bip,  cip, hc, 
                                                            K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, 
                                                            b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                        if ( gv.ABS_TYPE ==2 )
                            update_s_visc_vti_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pr, pp, pq,
                                                            pc55ipjpu, pc13u, pc11u,  pc33u, pc55ipjpd,  pc13d, pc11d,  pc33d,
                                                            bip,  cip, absorb_coeff, hc, &gv );
                    }
#ifdef EBUG
		    debug_check_matrix(psxx, nt, gv.NX, gv.NY, 666, 0, "psxx");
		    debug_check_matrix(psyy, nt, gv.NX, gv.NY, 666, 0, "psyy");
		    debug_check_matrix(psxy, nt, gv.NX, gv.NY, 666, 0, "psxy");
#endif
                    break;
                        
                        
                    case 7: /* elastic TTI */
                        update_s_elastic_TTI_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy,
                                                       pvyx, pvxy, psxx, psyy, psxy,
                                                       pc11, pc55ipjp, pc13, pc33, pc15, pc35, pc15ipjp, pc35ipjp, hc, &gv );
                       
                     if ( gv.FW ) {
                         if ( gv.ABS_TYPE ==1 )
                       update_s_elastic_TTI_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy,
                                                   pvyx, pvxy, psxx, psyy, psxy,
                                                   pc11, pc55ipjp, pc13, pc33, pc15, pc35, pc15ipjp, pc35ipjp, hc,
                                               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                    

                             
                         if ( gv.ABS_TYPE ==2 )
                             update_s_elastic_tti_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy,
                                                       pvyx, pvxy, psxx, psyy, psxy,
                                                       pc11, pc55ipjp, pc13, pc33, pc15, pc35, pc15ipjp, pc35ipjp, absorb_coeff, hc );
                     }

#ifdef EBUG
		     debug_check_matrix(psxx, nt, gv.NX, gv.NY, 777, 0, "psxx");
		     debug_check_matrix(psyy, nt, gv.NX, gv.NY, 777, 0, "psyy");
		     debug_check_matrix(psxy, nt, gv.NX, gv.NY, 777, 0, "psxy");
#endif
                     break;
                        
                    case 8: /* viscoelastic TTI */

                        update_s_visc_TTI_interior ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy, pvyx, pvxy, psxx, psyy, psxy, pr, pp, pq,
                                                                         pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                                                         pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                                                            bip,  cip, hc, &gv);
                   if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_visc_TTI_PML ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy, pvyx, pvxy, psxx, psyy, psxy, pr, pp, pq,
                                                   pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                                   pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                                   bip,  cip, hc, K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y,
                                                   K_y_half, a_y_half,b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx, &gv );
                        if ( gv.ABS_TYPE ==2 )
                            update_s_visc_tti_abs ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvxx, pvyy, pvyx, pvxy, psxx, psyy, psxy, pr, pp, pq,
                                                   pc11u, pc33u,  pc13u, pc55u, pc15u, pc35u, pc55ipjpu, pc15ipjpu, pc35ipjpu,
                                                   pc11d, pc33d,  pc13d, pc55d, pc15d, pc35d, pc55ipjpd, pc15ipjpd, pc35ipjpd,
                                                            bip,  cip, absorb_coeff, hc, &gv );
                    }
#ifdef EBUG
		    debug_check_matrix(psxx, nt, gv.NX, gv.NY, 888, 0, "psxx");
		    debug_check_matrix(psyy, nt, gv.NX, gv.NY, 888, 0, "psyy");
		    debug_check_matrix(psxy, nt, gv.NX, gv.NY, 888, 0, "psxy");
#endif
                    break;
                }
            }
            
            if(gv.FDORDER_TIME==4){
                if ( gv.L ) { /* viscoelastic */
                    /* Not supported right now */
                    update_s_visc_interior_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, pr, pp, pq, fipjp,f, g, bip, bjm, cip, cjm, d, e, dip, hc ,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4, &gv);
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 ) {
                            update_s_visc_PML_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, hc, pr, pp, pq, fipjp,
                                                 f, g, bip, bjm, cip, cjm, d, e, dip,K_x, a_x, b_x, K_x_half, a_x_half, b_x_half,
                                                 K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx ,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4, &gv);
                        }
                        if ( gv.ABS_TYPE !=1 ) {
                            update_s_visc_abs_4 ( 1, gv.NX, 1, gv.NY, gx,gy, nt, pvx, pvy, psxx, psyy, psxy, pr,
                                                 pp, pq, ppi, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
                                                 absorb_coeff,hc ,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4, &gv);
                        }
                    }
                    /* Shift memory variables one time step back */
                    shift_r1=pp_4;pp_4=pp_3;pp_3=pp_2;pp_2=pp;pp=shift_r1;
                    shift_r2=pr_4;pr_4=pr_3;pr_3=pr_2;pr_2=pr;pr=shift_r2;
                    shift_r3=pq_4;pq_4=pq_3;pq_3=pq_2;pq_2=pq;pq=shift_r3;
                } else { /* elastic */
                    update_s_elastic_interior_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, ppi, pu, puipjp, hc,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, &gv);
                    
                    if ( gv.FW ) {
                        if ( gv.ABS_TYPE ==1 )
                            update_s_elastic_PML_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy, ppi, pu, puipjp, hc,
                                                  K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, &gv);
                        if ( gv.ABS_TYPE !=1 )
                            update_s_elastic_abs_4 ( 1, gv.NX, 1, gv.NY, gx, gy, nt, pvx, pvy, psxx, psyy, psxy,
                                                  ppi, pu, puipjp, absorb_coeff, hc,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4, &gv);
                    }
                }
                /* Shift spartial derivations from the velocity one time step back */
                shift_v1=vxx_4;vxx_4=vxx_3;vxx_3=vxx_2;vxx_2=vxx_1;vxx_1=shift_v1;
                shift_v2=vyy_4;vyy_4=vyy_3;vyy_3=vyy_2;vyy_2=vyy_1;vyy_1=shift_v2;
                shift_v3=vxy_4;vxy_4=vxy_3;vxy_3=vxy_2;vxy_2=vxy_1;vxy_1=shift_v3;
                shift_v4=vyx_4;vyx_4=vyx_3;vyx_3=vyx_2;vyx_2=vyx_1;vyx_1=shift_v4;
            }
            
            /* explosive source */
            if ( gv.SOURCE_TYPE == 1 )
                psource ( nt, psxx, psyy, srcpos_loc, signals, nsrc_loc, &gv );
            
            if ( ( gv.FREE_SURF ) && ( gv.POS[2] == 0 ) ) {
                if ( gv.L ) /* viscoelastic */
                    surface ( 1, pvx, pvy, psxx, psyy, psxy, pp, pq, ppi, pu,
                             ptaup, ptaus, etajm, peta, hc, K_x, a_x, b_x, psi_vxx, &gv );
                else
                /* elastic */
                    surface_elastic ( 1, gx, pvx, pvy, psxx, psyy, psxy, ppi, pu, hc, K_x, a_x, b_x, psi_vxxs, &gv );
            }
            
            if ( ( gv.MPID == 0 )
                && ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
                time6 = MPI_Wtime();
                time_av_s_update += ( time6 - time5 );
                log_debug("Starting stress exchange between PEs...\n");
            }
            /*---------------------------------------------------------------*/
            /* -------- stress exchange between PEs --------*/
            /*---------------------------------------------------------------*/
             
            /*if ( RSG ) {
             exchange_s_rsg ( psxx, psyy, psxy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, &gv );
             } else {*/
            exchange_s (nd, psxx, psyy, psxy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec, &gv );
            
            if ( ( gv.MPID == 0 )	&& ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
                
                time7 = MPI_Wtime();
                time_av_s_exchange += ( time7 - time6 );
                log_debug("Finished stress exchange between PEs (real time: %4.3fs).\n", time7 - time6);
            }
            
            /* store amplitudes at receivers in section-arrays */
            if ( ( gv.SEISMO ) && ( nt == lsamp ) && ( nt < gv.NT ) ) {
                
                seismo_ssg ( lsamp, ntr, recpos_loc, sectionvx, sectionvy, sectionp, sectioncurl, sectiondiv, pvx, pvy, psxx, psyy, ppi, pu, hc, &gv );
                lsamp += gv.NDT;
            }
            
            /* WRITE SNAPSHOTS TO DISK */
            if ( ( gv.SNAP ) && ( nt == lsnap ) && ( nt <= gv.TSNAP2 / gv.DT ) ) {
                
                snap (nt, ++nsnap, pvx, pvy, psxx, psyy, pu, ppi, hc, &gv );
                lsnap = lsnap + iround ( gv.TSNAPINC/gv.DT );
            }
            
            if ( ( gv.MPID == 0 )	&& ( ( nt + ( gv.OUTNTIMESTEPINFO - 1 ) ) % gv.OUTNTIMESTEPINFO ) == 0 ) {
                
                time8 = MPI_Wtime();
                time_av_timestep += ( time8 - time3 );
                log_info("Total real time for timestep %d: %4.3fs.\n", nt, time8 - time3);
            }
        }
        /*---------------------------------------------------------------*/
        /*--------------------  End  of loop over timesteps ----------*/
        /*---------------------------------------------------------------*/
        
        log_infoc(0, "Finished time stepping.\n");
         
        /* write seismograms to file(s) */
        if ( gv.SEISMO ) {
            
            /* saves seismograms portion of each PE individually to file */
            //if (ntr> 0) saveseis(sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos_current,ishot,ns);
            
            /* merge of seismogram data from all PE and output data collectively */
            switch ( gv.SEISMO ) {
                case 1 : /* particle velocities only */
                    catseis ( sectionvx, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,1, &gv );
                    catseis ( sectionvy, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,2, &gv );
                    
                    break;
                case 2 : /* pressure only */
                    catseis ( sectionp, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,4, &gv );
                    
                    break;
                case 3 : /* curl and div only */
                    catseis ( sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,5, &gv );
                    catseis ( sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,6, &gv );
                    
                    break;
                case 4 : /* everything */
                    catseis ( sectionvx, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,1, &gv );
                    catseis ( sectionvy, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,2, &gv );
                    catseis ( sectionp, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,4, &gv );
                    catseis ( sectiondiv, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,5, &gv );
                    catseis ( sectioncurl, seismo_fulldata, recswitch, ntr_glob,ns );
                    if ( gv.MPID==0 ) saveseis_glob (seismo_fulldata,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,6, &gv );
                    
                    break;
                default :
                    break;
                    
            }
	} 
    } /* end of loop over shots */
    
    /* deallocation of memory */
    free_matrix ( psxx, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( psxy, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( psyy, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( pvx, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( pvy, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    
    if (gv.FDORDER_TIME==4) {
        free_matrix (vxx_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxx_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxx_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxx_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        free_matrix (vyy_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyy_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyy_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyy_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        free_matrix (vxy_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxy_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxy_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vxy_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        free_matrix (vyx_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyx_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyx_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (vyx_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        free_matrix (svx_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svx_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svx_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svx_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        
        free_matrix (svy_1, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svy_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svy_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix (svy_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    }

   if (gv.WEQ>2){ 
    free_matrix ( prho, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( prip, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( prjp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( ppi, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( pu, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( puipjp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
    free_matrix ( absorb_coeff, 1, gv.NY, 1, gv.NX );
}
    
    free_ivector ( gx,1,4 );
    free_ivector ( gy,1,4 );
    
    if ( gv.L ) {
        free_f3tensor ( pr, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pq, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
    }


if ( gv.WEQ==4 ) { /*viscoelastic wave equation */
        free_matrix ( ptaus, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix ( ptausipjp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix ( ptaup, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_vector ( peta, 1, gv.L );
        free_vector ( etaip, 1, gv.L );
        free_vector ( etajm, 1, gv.L );
        free_vector ( bip, 1, gv.L );
        free_vector ( bjm, 1, gv.L );
        free_vector ( cip, 1, gv.L );
        free_vector ( cjm, 1, gv.L );
        free_matrix ( f, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix ( g, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_matrix ( fipjp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
        free_f3tensor ( dip, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( d, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( e, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
    }
    
   if ( gv.WEQ==6 ) { /*viscoelastic VTI wave equation */

            free_matrix ( pc33, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( pc13, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( ptau11, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( ptau33, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( ptau13, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( ptau55, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( ptau55ipjp, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

            free_f3tensor ( pc55ipjpd, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
            free_f3tensor ( pc13d, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
            free_f3tensor ( pc33d, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
            free_f3tensor ( pc11d, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );


            free_matrix ( pc55ipjpu, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( pc13u, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( pc11u, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );
            free_matrix ( pc33u, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd );

            free_vector ( peta, 1, gv.L );
            free_vector ( bip, 1, gv.L );
            free_vector ( cip, 1, gv.L );
}
    /*elastic TTI wave equation */
/*    if ( WEQ==7 ) {
        free_matrix ( pvxx, -nd + 1, NY + nd, -nd + 1, NX + nd );
        free_matrix ( pvyy, -nd + 1, NY + nd, -nd + 1, NX + nd );
        free_matrix ( pvxy, -nd + 1, NY + nd, -nd + 1, NX + nd );
        free_matrix ( pvyx, -nd + 1, NY + nd, -nd + 1, NX + nd );

    }*/

    if(gv.L>0 && gv.FDORDER_TIME==4){
        free_f3tensor ( pr_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pr_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pr_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        
        free_f3tensor ( pp_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pp_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pp_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        
        free_f3tensor ( pq_2, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pq_3, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
        free_f3tensor ( pq_4, -nd + 1, gv.NY + nd, -nd + 1, gv.NX + nd, 1, gv.L );
    }
    /*if ( RSG ) {
     free_matrix ( bufferlef_to_rig, 0, gv.NY + 1, 1, fdo3 );
     free_matrix ( bufferrig_to_lef, 0, gv.NY + 1, 1, fdo3 );
     free_matrix ( buffertop_to_bot, 1, gv.NX, 1, fdo3 );
     free_matrix ( bufferbot_to_top, 1, gv.NX, 1, fdo3 );
     } else {*/
    free_matrix ( bufferlef_to_rig, 1, gv.NY, 1, fdo3 );
    free_matrix ( bufferrig_to_lef, 1, gv.NY, 1, fdo3 );
    free_matrix ( buffertop_to_bot, 1, gv.NX, 1, fdo3 );
    free_matrix ( bufferbot_to_top, 1, gv.NX, 1, fdo3 );
    
    
    if ( nsrc_loc > 0 ) {
        if (signals) free_matrix ( signals, 1, nsrc_loc, 1, gv.NT );
        free_matrix ( srcpos_loc, 1, 12, 1, nsrc_loc );
    }
    
    if ( gv.ABS_TYPE==1 ) {
        
        free_vector ( d_x,1,2*gv.FW );
        free_vector ( K_x,1,2*gv.FW );
        free_vector ( alpha_prime_x,1,2*gv.FW );
        free_vector ( a_x,1,2*gv.FW );
        free_vector ( b_x,1,2*gv.FW );
        
        free_vector ( d_x_half,1,2*gv.FW );
        free_vector ( K_x_half,1,2*gv.FW );
        free_vector ( alpha_prime_x_half,1,2*gv.FW );
        free_vector ( a_x_half,1,2*gv.FW );
        free_vector ( b_x_half,1,2*gv.FW );
        
        free_vector ( d_y,1,2*gv.FW );
        free_vector ( K_y,1,2*gv.FW );
        free_vector ( alpha_prime_y,1,2*gv.FW );
        free_vector ( a_y,1,2*gv.FW );
        free_vector ( b_y,1,2*gv.FW );
        
        free_vector ( d_y_half,1,2*gv.FW );
        free_vector ( K_y_half,1,2*gv.FW );
        free_vector ( alpha_prime_y_half,1,2*gv.FW );
        free_vector ( a_y_half,1,2*gv.FW );
        free_vector ( b_y_half,1,2*gv.FW );
        
        free_matrix ( psi_sxx_x,1,gv.NY,1,2*gv.FW );
        free_matrix ( psi_syy_y,1,2*gv.FW,1,gv.NX );
        free_matrix ( psi_sxy_x,1,gv.NY,1,2*gv.FW );
        free_matrix ( psi_sxy_y,1,2*gv.FW,1,gv.NX );
        free_matrix ( psi_vxx,1,gv.NY,1,2*gv.FW );
        free_matrix ( psi_vyy,1,2*gv.FW,1,gv.NX );
        free_matrix ( psi_vxy,1,2*gv.FW,1,gv.NX );
        free_matrix ( psi_vyx,1,gv.NY,1,2*gv.FW );
        
    }
    
    
    
    if ( gv.SEISMO )
        free_imatrix ( recpos, 1, 3, 1, ntr_glob );
    
    /* free memory for global source positions */
    free_matrix ( srcpos, 1, 8, 1, nsrc );
    
    if ( ( ntr > 0 ) && ( gv.SEISMO ) ) {
        
        free_matrix ( seismo_fulldata,1,ntr_glob,1,ns );
        free_imatrix ( recpos_loc, 1, 3, 1, ntr );
        switch ( gv.SEISMO ) {
            case 1: /* particle velocities only */
                free_matrix ( sectionvx, 1, ntr, 1, ns );
                free_matrix ( sectionvy, 1, ntr, 1, ns );
                break;
            case 2: /* pressure only */
                free_matrix ( sectionp, 1, ntr, 1, ns );
                break;
            case 3: /* curl and div only */
                free_matrix ( sectioncurl, 1, ntr, 1, ns );
                free_matrix ( sectiondiv, 1, ntr, 1, ns );
                break;
            case 4: /* everything */
                free_matrix ( sectionvx, 1, ntr, 1, ns );
                free_matrix ( sectionvy, 1, ntr, 1, ns );
                free_matrix ( sectionp, 1, ntr, 1, ns );
                free_matrix ( sectioncurl, 1, ntr, 1, ns );
                free_matrix ( sectiondiv, 1, ntr, 1, ns );
                break;
        }    
    }
    
    /* de-allocate buffer for messages */
    MPI_Buffer_detach(buff_addr, &buffsize);

    MPI_Barrier (MPI_COMM_WORLD);

    if (gv.MPID == 0) {
      time_av_v_update = time_av_v_update / (double)gv.NT;
      time_av_s_update = time_av_s_update / (double)gv.NT;
      time_av_v_exchange = time_av_v_exchange / (double)gv.NT;
      time_av_s_exchange = time_av_s_exchange / (double)gv.NT;
      time_av_timestep = time_av_timestep / (double)gv.NT;
      log_info("Average times for\n");
      log_info("  velocity update: .. %5.6fs.\n", time_av_v_update);
      log_info("  stress update: .... %5.6fs.\n", time_av_s_update);
      log_info("  velocity exchange:  %5.6fs.\n", time_av_v_exchange);
      log_info("  stress exchange: .. %5.6fs.\n", time_av_s_exchange);
      log_info("  time step: ........ %5.6fs.\n", time_av_timestep);
      cpu_time = clock()-cpu_time1;
      log_info("CPU time of program per PE: %lis.\n",cpu_time/CLOCKS_PER_SEC);
      time8 = MPI_Wtime();
      log_info("Total real time of program: %4.3fs.\n", time8-time1);
    }
    
    /* finalize logging */
    log_finalize();

    if (log_fp) fclose(log_fp);
    
    /* finalize MPI */
    MPI_Finalize();

    return EXIT_SUCCESS;
}
