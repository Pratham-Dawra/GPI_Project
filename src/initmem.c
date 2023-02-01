
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

/* -------------------------------------------------------------
 * This is function initmem.
 * Initialising memory for variables.
 *
 * -------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void initmem(MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{

    int nseismograms = 0;

    float memdyn, memmodel, memseismograms, membuffer, memtotal, memcpml = 0.0;
    float fac1, fac2, memadd, memadd_L;

    /* number of seismogram sections which have to be stored in core memory */
    /* allocation of memory for seismogramm merge */
    switch (gv->SEISMO) {
      case 1:                  /* particle velocities only */
          nseismograms = 3;
          break;
      case 2:                  /* pressure only */
          nseismograms = 1;
          break;
      case 3:                  /* curl and div only */
          nseismograms = 2;
          break;
      case 4:                  /* everything */
          nseismograms = 6;
          break;
      default:
          nseismograms = 1;
          break;
    }

    /*estimate memory requirements for dynamic, static and buffer arrays */
    gv->ND = gv->FDORDER / 2;
    gv->FDO3 = 2 * gv->ND;      // WHY ????

    fac1 = (gv->NX + gv->FDO3) * (gv->NY + gv->FDO3);
    fac2 = sizeof(float) * pow(2.0, -20.0);
    memadd = 0.0;
    memadd_L = 0.0;
    if (gv->FDORDER_TIME == 4) {
        memadd = 24.0;
        memadd_L = 9;
    }
    if (gv->L) {
        memdyn = (memadd + 5.0 + (3.0 + memadd_L) * (float)gv->L) * fac1 * fac2;
        memmodel = (12.0 + 3.0 * (float)gv->L) * fac1 * fac2 + gv->NX * gv->NY * fac2;
    } else {
        memdyn = (memadd + 5.0) * fac1 * fac2;
        memmodel = 6.0 * fac1 * fac2 + gv->NX * gv->NY * fac2;
    }

    memseismograms = nseismograms * gv->NTR * gv->NS * fac2;
    membuffer = 2.0 * gv->FDO3 * (gv->NY + gv->NX) * fac2;
    gv->BUFFSIZE = 2.0 * 2.0 * gv->FDO3 * (gv->NX + gv->NY) * sizeof(MPI_FLOAT);
    if (gv->ABS_TYPE == 1)
        memcpml = 2.0 * gv->FW * 4.0 * (gv->NY + gv->NX) * fac2 + 20.0 * 2.0 * gv->FW * fac2;
    memtotal = memdyn + memmodel + memseismograms + membuffer + memcpml + (gv->BUFFSIZE * pow(2.0, -20.0));

    if (gv->MPID == 0) {
        log_info("Size of local grids: NX=%d, NY=%d\n", gv->NX, gv->NY);
        log_info("Each process is now trying to allocate memory for:\n");
        log_info("  Dynamic variables: ............. %6.2f MB\n", memdyn);
        log_info("  Static variables: .............. %6.2f MB\n", memmodel);
        log_info("  Seismograms: ................... %6.2f MB\n", memseismograms);
        log_info("  Buffer arrays for grid exchange: %6.2f MB\n", membuffer);
        log_info("  Network buffer for MPI_Bsend: .. %6.2f MB\n", gv->BUFFSIZE * pow(2.0, -20.0));
        if (gv->ABS_TYPE == 1)
            log_info("  CPML variables: ................ %6.2f MB\n", memcpml);
        log_info("------------------------------------------------\n");
        log_info("Total memory required: ........... %6.2f MB.\n", memtotal);
    }

    /* allocate buffer for buffering messages */
    gv->BUFF_ADDR = malloc(gv->BUFFSIZE);
    if (!gv->BUFF_ADDR)
        log_fatal("Allocation failure for MPI_Bsend buffer!\n");
    MPI_Buffer_attach(gv->BUFF_ADDR, gv->BUFFSIZE);

    /* ------------ memory allocation for arrays ------------- */
    /* subgrid arrays */
    gv->GY = ivector(1, 4);
    gv->GX = ivector(1, 4);

    /* Memory allocation for wavefields */
    initmem_wavefield(mpw, gv);

    /* Memory allocation for models */
    initmem_model(mpm, gv);

    /* allocate buffer for seismogram output, merged seismogram section of all PEs */
    if (gv->SEISMO)
        gv->SEISMO_FULLDATA = matrix(1, gv->NTRG, 1, gv->NS);

    if (gv->NTR > 0) {
        switch (gv->SEISMO) {
          case 1:              /* particle velocities only */
              gv->SECTIONVX = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONVY = matrix(1, gv->NTR, 1, gv->NS);
              break;
          case 2:              /* pressure only */
              gv->SECTIONP = matrix(1, gv->NTR, 1, gv->NS);
              break;
          case 3:              /* curl and div only */
              gv->SECTIONCURL = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONDIV = matrix(1, gv->NTR, 1, gv->NS);
              break;
          case 4:              /* everything */
              gv->SECTIONVX = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONVY = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONCURL = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONDIV = matrix(1, gv->NTR, 1, gv->NS);
              gv->SECTIONP = matrix(1, gv->NTR, 1, gv->NS);
              break;
        }
    }

    log_debug("Memory allocation for PE %d was successful.\n", gv->MPID);

}
