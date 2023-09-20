
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * inversion process -
 *----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void inversion(int iter, int ishot, AcqVar *acq, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    int h, i, j;
    int itestshot, swstestshot;
    float **sectionread = NULL, **fulldata = NULL;

    sectionread = matrix(1, gv->NTRG, 1, gv->NS);
    if (vinv->WRITE_DIFF || vinv->WRITE_FILTERED_DATA)
        fulldata = matrix(1, gv->NTRG, 1, gv->NS);

    /*-----------------------------------*/

    /*------- Calculate residuals -------*/

    /*-----------------------------------*/

    if (gv->MPID == 0) {
        log_info("---------------------------\n");
        log_info("--- Calculate residuals ---\n");
        log_info("---------------------------\n");
    }

    itestshot = vinv->TESTSHOT_START;

    log_info("MPID, NTR: %d %d\n", gv->MPID, gv->NTR);
    if (gv->NTR > 0) {
        log_info("MPID, NTR: %d %d\n", gv->MPID, gv->NTR);

        /* calculate L2-Norm and energy ? */
        if ((ishot == itestshot) && (ishot <= vinv->TESTSHOT_END)) {
            swstestshot = 1;
        }
    }

    /* --------------------------------------------- */
    /* ----- read seismic data from SU file vx ----- */
    /* --------------------------------------------- */
    /* Output test */
    float **srcpos = matrix(1, NSPAR, 1, 1);
    int **recpos = imatrix(1, 3, 1, gv->NTRG);

    if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {   /* if ADJOINT_TYPE */
        inseis(ishot, 1, sectionread, gv, vinv);
        /* Output test */
        outseis_glob(fopen("Test_seisread_prefilt_vx.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        //log_info("TIME_FILT; %d, FREQ_NR: %d\n", vinv->TIME_FILT, vinv->FREQ_NR);
        if ((vinv->TIME_FILT == 1) || (vinv->TIME_FILT == 2)) {
            timedomain_filt(sectionread, 1, gv->NTRG, gv, vinv);
            timedomain_filt(minv->sectionvxcalc, 1, gv->NTR, gv, vinv);
        }
        outseis_glob(fopen("Test_seisread_postfilt_vx.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        h = 1;
        for (i = 1; i <= gv->NTR; i++) {
            for (j = 1; j <= gv->NS; j++) {
                minv->sectionvxdata[h][j] = sectionread[acq->recpos_loc[3][i]][j];
            }
            h++;
        }
        calc_res(minv->sectionvxdata, minv->sectionvxcalc, minv->sectionvxdiff, minv->sectionvxdiffold, 1, swstestshot,
                 ishot, iter, acq, gv, vinv);
        /*if (swstestshot == 1) {
         * energy =
         * calc_energy(sectionvxdata, ntr, ns, energy, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * }
         * energy_all_shots =
         * calc_energy(sectionvxdata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * L2_all_shots =
         * calc_misfit(sectionvxdata, sectionvx, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
         * ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos); */
        /*fprintf(FP,"Energy vxdata for PE %d:   %f\n\n", MYID,energy); */
    }

    /* end ADJOINT_TYPE */

    /* --------------------------------- */
    /* read seismic data from SU file vy */
    /* --------------------------------- */
    if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {   /* if ADJOINT_TYPE */
        inseis(ishot, 2, sectionread, gv, vinv);
        outseis_glob(fopen("Test_seisread_prefilt_vy.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        if ((vinv->TIME_FILT == 1) || (vinv->TIME_FILT == 2)) {
            timedomain_filt(sectionread, 1, gv->NTRG, gv, vinv);
            timedomain_filt(minv->sectionvycalc, 1, gv->NTR, gv, vinv);
        }
        outseis_glob(fopen("Test_seisread_postfilt_vy.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        h = 1;
        for (i = 1; i <= gv->NTR; i++) {
            for (j = 1; j <= gv->NS; j++) {
                minv->sectionvydata[h][j] = sectionread[acq->recpos_loc[3][i]][j];
            }
            h++;
        }
        calc_res(minv->sectionvydata, minv->sectionvycalc, minv->sectionvydiff, minv->sectionvydiffold, 1, swstestshot,
                 ishot, iter, acq, gv, vinv);

        /*if (swstestshot == 1) {
         * energy =
         * calc_energy(sectionvydata, ntr, ns, energy, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * }
         * energy_all_shots =
         * calc_energy(sectionvydata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * L2_all_shots =
         * calc_misfit(sectionvydata, sectionvy, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
         * ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos); */
        /*fprintf(FP,"Energy vydata for PE %d:   %f\n\n", MYID,energy); */
    }

    /* end ADJOINT_TYPE */

    /* --------------------------------- */
    /* read seismic data from SU file p */
    /* --------------------------------- */
    if (vinv->ADJOINT_TYPE == 4) {  /* if ADJOINT_TYPE */
        inseis(ishot, 9, sectionread, gv, vinv);
        outseis_glob(fopen("Test_seisread_prefilt_p.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        if ((vinv->TIME_FILT == 1) || (vinv->TIME_FILT == 2)) {
            timedomain_filt(sectionread, 1, gv->NTRG, gv, vinv);
            timedomain_filt(minv->sectionpcalc, 1, gv->NTR, gv, vinv);
        }
        outseis_glob(fopen("Test_seisread_postfilt_p.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS,
                     gv->SEIS_FORMAT, ishot, 1, gv);
        h = 1;
        for (i = 1; i <= gv->NTR; i++) {
            for (j = 1; j <= gv->NS; j++) {
                minv->sectionpdata[h][j] = sectionread[acq->recpos_loc[3][i]][j];
            }
            h++;
        }
        calc_res(minv->sectionpdata, minv->sectionpcalc, minv->sectionpdiff, minv->sectionpdiffold, 1, swstestshot,
                 ishot, iter, acq, gv, vinv);
        /*if (swstestshot == 1) {
         * energy =
         * calc_energy(sectionpdata, ntr, ns, energy, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * }
         * energy_all_shots =
         * calc_energy(sectionpdata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
         * nsrc_glob, ishot, iter, srcpos, recpos);
         * L2_all_shots =
         * calc_misfit(sectionpdata, sectionp, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
         * ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos); */
    }

    /* end ADJOINT_TYPE */
    /* --------------------------------- */

    // Tracekill
    /*if (TRKILL) {
     * count_killed_traces(ntr, swstestshot, ntr_glob, recpos_loc, nsrc_glob, ishot, ptr_killed_traces, ptr_killed_traces_testshots, srcpos, recpos);
     * } */

    if ((ishot == itestshot) && (ishot <= vinv->TESTSHOT_END)) {
        swstestshot = 0;
        itestshot += vinv->TESTSHOT_INCR;
    }


    /* ------------------------------------------- */
    /* write filtered data and differences to file */
    /* ------------------------------------------- */

    /* Write differences between measured and synthetic seismogramms (adjoint sources) to disk */
    if (vinv->WRITE_DIFF) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxdiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 11, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvydiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 12, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpdiff, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 14, gv);
        }
    }
    /* Write measured filtered seismogramms to disk */
    if (vinv->TIME_FILT && vinv->WRITE_FILTERED_DATA) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxdata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 21, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvydata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 22, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpdata, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 24, gv);
        }
    }
    /* Write synthetic filtered seismogramms to disk */
    if (vinv->TIME_FILT && vinv->WRITE_FILTERED_DATA == 2) {
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {
            catseis(minv->sectionvxcalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 31, gv);
        }
        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 2)) {
            catseis(minv->sectionvycalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 32, gv);
        }
        if (vinv->ADJOINT_TYPE == 4) {
            catseis(minv->sectionpcalc, fulldata, acq->recswitch, gv->NTRG, gv->NS);
            if (gv->MPID == 0)
                saveseis_glob(fulldata, acq->recpos, acq->srcpos, ishot, gv->NS, 34, gv);
        }
    }
    
    /* ------------------------------------------- */
    
    free_matrix(sectionread, 1, gv->NTRG, 1, gv->NS);
    if (vinv->WRITE_DIFF || vinv->WRITE_FILTERED_DATA)
        free_matrix(fulldata, 1, gv->NTRG, 1, gv->NS);

}
