
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

void inversion(int iter, int ishot, GlobVar *gv, GlobVarInv *vinv)
{
    int itestshot, swstestshot;
    float **sectionread;
    
    sectionread = matrix(1, gv->NTRG, 1, gv->NS+1);

    /*-----------------------------------*/
    /*------- Calculate residuals -------*/
    /*-----------------------------------*/

    if (gv->MPID == 0) {
        log_info("---------------------------\n");
        log_info("--- Calculate residuals ---\n");
        log_info("---------------------------\n");
    }

    itestshot = vinv->TESTSHOT_START;

    log_info("MPID, NTR: %d %d\n",gv->MPID, gv->NTR);
    if (gv->NTR > 0) {
        log_info("MPID, NTR: %d %d\n",gv->MPID, gv->NTR);

        /* calculate L2-Norm and energy ? */
        if ((ishot == itestshot) && (ishot <= vinv->TESTSHOT_END)) {
            swstestshot = 1;
        }
    }

    /* --------------------------------------------- */
    /* ----- read seismic data from SU file vx ----- */
    /* --------------------------------------------- */

        if ((vinv->ADJOINT_TYPE == 1) || (vinv->ADJOINT_TYPE == 3)) {   /* if ADJOINT_TYPE */
            inseis(iter, ishot, 1, sectionread, gv, vinv);
            float **srcpos = matrix(1, NSPAR, 1, 1);
            int **recpos = imatrix(1, 3, 1, gv->NTRG);
            outseis_glob(fopen("Test_seisread_3.su", "w"), sectionread, recpos, gv->NTRG, srcpos, gv->NS, gv->SEIS_FORMAT, ishot, 1, gv);
            /*if ((TIME_FILT == 1) || (TIME_FILT == 2)) {
                timedomain_filt(sectionread, F_LOW_PASS, ORDER, ntr_glob, ns, 1);
            }*/
            /*h = 1;
            for (i = 1; i <= ntr; i++) {
                for (j = 1; j <= ns; j++) {
                    sectionvxdata[h][j] = sectionread[recpos_loc[3][i]][j];
                }
                h++;
            }
            L2 = calc_res(sectionvxdata, sectionvx, sectionvxdiff, sectionvxdiffold, ntr,
                          ns, LNORM, L2, 0, 1, swstestshot, ntr_glob, recpos_loc, nsrc_glob,
                          ishot, iter, srcpos, recpos);
            if (swstestshot == 1) {
                energy =
                    calc_energy(sectionvxdata, ntr, ns, energy, ntr_glob, recpos_loc,
                                nsrc_glob, ishot, iter, srcpos, recpos);
            }
            energy_all_shots =
                calc_energy(sectionvxdata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
                            nsrc_glob, ishot, iter, srcpos, recpos);
            L2_all_shots =
                calc_misfit(sectionvxdata, sectionvx, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
                            ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos);*/
            /*fprintf(FP,"Energy vxdata for PE %d:   %f\n\n", MYID,energy); */
        }

        /* end ADJOINT_TYPE */
        /* --------------------------------- */
        /* read seismic data from SU file vy */
        /* --------------------------------- */
        //if ((ADJOINT_TYPE == 1) || (ADJOINT_TYPE == 2)) {   /* if ADJOINT_TYPE */
        /*    inseis(fprec, ishot, sectionread, ntr_glob, ns, 2, iter);
            if ((TIME_FILT == 1) || (TIME_FILT == 2)) {
                timedomain_filt(sectionread, F_LOW_PASS, ORDER, ntr_glob, ns, 1);
            }
            h = 1;
            for (i = 1; i <= ntr; i++) {
                for (j = 1; j <= ns; j++) {
                    sectionvydata[h][j] = sectionread[recpos_loc[3][i]][j];
                }
                h++;
            }
            L2 = calc_res(sectionvydata, sectionvy, sectionvydiff, sectionvydiffold, ntr,
                          ns, LNORM, L2, 0, 1, swstestshot, ntr_glob, recpos_loc, nsrc_glob,
                          ishot, iter, srcpos, recpos);
            if (swstestshot == 1) {
                energy =
                    calc_energy(sectionvydata, ntr, ns, energy, ntr_glob, recpos_loc,
                                nsrc_glob, ishot, iter, srcpos, recpos);
            }
            energy_all_shots =
                calc_energy(sectionvydata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
                            nsrc_glob, ishot, iter, srcpos, recpos);
            L2_all_shots =
                calc_misfit(sectionvydata, sectionvy, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
                            ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos);*/
            /*fprintf(FP,"Energy vydata for PE %d:   %f\n\n", MYID,energy); */
        //}

        /* end ADJOINT_TYPE */
        /* --------------------------------- */
        /* read seismic data from SU file p */
        /* --------------------------------- */
        //if (ADJOINT_TYPE == 4) {    /* if ADJOINT_TYPE */
        /*    inseis(fprec, ishot, sectionread, ntr_glob, ns, 9, iter);
            if ((TIME_FILT == 1) || (TIME_FILT == 2)) {
                timedomain_filt(sectionread, F_LOW_PASS, ORDER, ntr_glob, ns, 1);
            }
            h = 1;
            for (i = 1; i <= ntr; i++) {
                for (j = 1; j <= ns; j++) {
                    sectionpdata[h][j] = sectionread[recpos_loc[3][i]][j];
                }
                h++;
            }
            L2 = calc_res(sectionpdata, sectionp, sectionpdiff, sectionpdiffold, ntr, ns,
                          LNORM, L2, 0, 1, swstestshot, ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos);
            if (swstestshot == 1) {
                energy =
                    calc_energy(sectionpdata, ntr, ns, energy, ntr_glob, recpos_loc,
                                nsrc_glob, ishot, iter, srcpos, recpos);
            }
            energy_all_shots =
                calc_energy(sectionpdata, ntr, ns, energy_all_shots, ntr_glob, recpos_loc,
                            nsrc_glob, ishot, iter, srcpos, recpos);
            L2_all_shots =
                calc_misfit(sectionpdata, sectionp, ntr, ns, LNORM, L2_all_shots, 0, 1, 1,
                            ntr_glob, recpos_loc, nsrc_glob, ishot, iter, srcpos, recpos);*/
        //}                       /* end ADJOINT_TYPE */
}
