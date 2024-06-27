
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

void inversion(int iter, int ishot, int snapcheck, float *hc, AcqVar *acq, MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVar *gv, GlobVarInv *vinv, Perform *perf)
{
    int h, i, j;
    int itestshot, swstestshot;
    float **sectionread = NULL;
    float **srcpos_loc_back = NULL;

    sectionread = matrix(1, gv->NTRG, 1, gv->NS);

    /*-----------------------------------*/

    /*------- Calculate residuals -------*/

    /*-----------------------------------*/

    if (gv->MPID == 0) {
        log_info("---------------------------\n");
        log_info("--- Calculate residuals ---\n");
        log_info("---------------------------\n");
    }

    itestshot = vinv->TESTSHOT_START;

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
    if (vinv->WRITE_DIFF || (vinv->TIME_FILT && vinv->WRITE_FILTERED_DATA))
        saveseis_fwi(ishot, acq, minv, gv, vinv);

    /* ------------------------------------------- */
    /*          starting backpropagation           */
    /* ------------------------------------------- */
    if (gv->MPID == 0) {
        log_info("============================================================\n");
        log_info("**********  Starting simulation (backward model)  **********\n");
        log_info("============================================================\n");
    }

    /* determine source position out of reciever position on grid */

    /* Distribute multiple source positions on subdomains */
    /* define source positions at the receivers */
    srcpos_loc_back = matrix(1, 6, 1, gv->NTR);
    for (i = 1; i <= gv->NTR; i++) {
        srcpos_loc_back[1][i] = (acq->recpos_loc[1][i]);
        srcpos_loc_back[2][i] = (acq->recpos_loc[2][i]);
    }
    //ntr1 = gv->NTR;

    /* initialize wavefield with zero   */
    zero_wavefield(iter, mpw, minv, gv, vinv);

    /*--------------------------------------------------------------------------------*/
    /*---------------------- Start loop over timesteps (backpropagation) -------------*/
    time_loop(iter, ishot, snapcheck, hc, srcpos_loc_back, minv->sectionvxdiff, minv->sectionvydiff,
              gv->NTR, 1, acq, mpm, mpw, minv, gv, vinv, perf);





    free_matrix(sectionread, 1, gv->NTRG, 1, gv->NS);

}
