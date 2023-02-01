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
--------------------------------------------------------------------------
 * memory allocation
--------------------------------------------------------------------------*/

#include "fd.h"

void **alloc_mem(GlobVar *gv){

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
    
    if (MYID == 0) {
        fprintf(gv.FP, "\n **Message from main (printed by PE %d):\n", MYID);
        fprintf(gv.FP, " Size of local grids: NX=%d \t NY=%d\n", gv.NX, gv.NY);
        fprintf(gv.FP, " Each process is now trying to allocate memory for:\n");
        fprintf(gv.FP, " Dynamic variables: \t\t %6.2f MB\n", memdyn);
        fprintf(gv.FP, " Static variables: \t\t %6.2f MB\n", memmodel);
        fprintf(gv.FP, " Seismograms: \t\t\t %6.2f MB\n", memseismograms);
        fprintf(gv.FP, " Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
        fprintf(gv.FP, " Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
        if (gv.ABS_TYPE == 1) fprintf(gv.FP, " CPML variables: \t\t %6.2f MB\n", memcpml);
        fprintf(gv.FP, " ------------------------------------------------ \n");
        fprintf(gv.FP, " Total memory required: \t %6.2f MB.\n\n", memtotal);
	fflush(gv.FP);
    }

    /* allocate buffer for buffering messages */
    buff_addr = malloc ( buffsize );
    if ( !buff_addr )
        declare_error ( "allocation failure for buffer for MPI_Bsend !" );
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

    }
    
    /* memory allocation for buffer arrays in which the wavefield
     information to be exchanged between neighboring PEs is stored */
    
    /*if ( RSG ) {
     // in the RSG case fdo3 is always 4
     bufferlef_to_rig = matrix ( 0, gv.NY + 1, 1, fdo3 );
     bufferrig_to_lef = matrix ( 0, gv.NY + 1, 1, fdo3 );
     buffertop_to_bot = matrix ( 1, gv.NX, 1, fdo3 );
     bufferbot_to_top = matrix ( 1, gv.NX, 1, fdo3 );
     
     } else {*/
    
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
    
    fprintf ( gv.FP, " ... memory allocation for PE %d was successfull.\n\n", MYID );
