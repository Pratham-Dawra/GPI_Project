
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

/*------------------------------------------------------------------------
 *   Calculate Data Residuals
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "logging.h"

void calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int sws,
              int swstestshot, int ishot, int iter, AcqVar *acq, GlobVar *gv, GlobVarInv *vinv)
{

    /* declaration of variables */
    float RMS, signL1 = 0.0;
    int h, i, j, invtime, umax = 0;
    float abs_section = 0.0, abs_sectiondata = 0.0, sectiondata_mult_section = 0.0;
    float intseis_s, intseis_sd;
    float *picked_times = NULL;
    float **intseis_section = NULL, **intseis_sectiondata = NULL;
    float **intseis_sectiondata_envelope = NULL, **intseis_section_envelope = NULL, **intseis_section_hilbert =
        NULL, **dummy_1 = NULL, **dummy_2 = NULL;

    if (vinv->LNORM == 8) {
        intseis_section_envelope = matrix(1, gv->NTR, 1, gv->NS);
        intseis_sectiondata_envelope = matrix(1, gv->NTR, 1, gv->NS);
        intseis_section_hilbert = matrix(1, gv->NTR, 1, gv->NS);
        dummy_1 = matrix(1, gv->NTR, 1, gv->NS);
        dummy_2 = matrix(1, gv->NTR, 1, gv->NS);
    }

    intseis_section = matrix(1, gv->NTR, 1, gv->NS);    /* declaration of variables for integration */
    intseis_sectiondata = matrix(1, gv->NTR, 1, gv->NS);
    if (vinv->TIMEWIN)
        picked_times = vector(1, gv->NTR);  /* declaration of variables for TIMEWIN */
    /* sectiondiff will be set to zero */
    //umax=gv->NTR*gv->NS;
    //zero(&sectiondiff[1][1],umax);
    umax = gv->NTR * gv->NS * sizeof(float);
    //bzero(&(sectiondiff[1][1]), umax);
    /* --------------------------------- TRACE KILLING --------------------------------- */
    /* declaration of variables for trace killing */
    int **kill_tmp = NULL, *kill_vector = NULL;
    char trace_kill_file[STRING_SIZE * 2];
    FILE *ftracekill;

    /* Start Tracekill */
    if (vinv->TRKILL) {
        kill_tmp = imatrix(1, gv->NTRG, 1, acq->nsrc);
        kill_vector = ivector(1, gv->NTR);

        /* clear kill table */
        for (i = 1; i <= acq->nsrc; i++) {
            for (j = 1; j <= gv->NTRG; j++) {
                kill_tmp[j][i] = 0;
            }
        }

        /* Use ONLY offset based TraceKill */
        if (vinv->TRKILL_OFFSET == 1) {
            /* Generate TraceKill file on the fly based on the offset from the source */
            //create_trkill_table(kill_tmp,gv->NTRG,acq->recpos,acq->nsrc,acq->srcpos,ishot,vinv->TRKILL_OFFSET_LOWER,vinv->TRKILL_OFFSET_UPPER);
        } else {
            /* READ TraceKill file from disk */
            if (vinv->USE_WORKFLOW) {
                sprintf(trace_kill_file, "%s_%i.dat", vinv->TRKILL_FILE, vinv->WORKFLOW_STAGE);
                ftracekill = fopen(trace_kill_file, "r");
                if (ftracekill == NULL) {
                    /* If Workflow TraceKill file not found use File without workflow extensions */
                    sprintf(trace_kill_file, "%s.dat", vinv->TRKILL_FILE);
                    ftracekill = fopen(trace_kill_file, "r");
                    if (ftracekill == NULL) {
                        log_fatal(" Trace kill file could not be opened!");
                    }
                }
            } else {
                sprintf(trace_kill_file, "%s.dat", vinv->TRKILL_FILE);
                ftracekill = fopen(trace_kill_file, "r");
                if (ftracekill == NULL) {
                    log_fatal(" Trace kill file could not be opened!");
                }
            }

            for (i = 1; i <= gv->NTRG; i++) {
                for (j = 1; j <= acq->nsrc; j++) {
                    fscanf(ftracekill, "%d", &kill_tmp[i][j]);
                    if (feof(ftracekill)) {
                        log_fatal("Error while reading TraceKill file. Check dimensions!\n");
                    }
                }
            }

            fclose(ftracekill);

            /* Use Tracekill FILE and add the offset based TraceKill */
            if (vinv->TRKILL_OFFSET == 2) {
                //create_trkill_table(kill_tmp,gv->NTRG,acq->recpos,acq->nsrc,acq->srcpos,-100,vinv->TRKILL_OFFSET_LOWER,vinv->TRKILL_OFFSET_UPPER);
            }
        }

        /* Generate local kill vector */
        h = 1;
        for (i = 1; i <= gv->NTR; i++) {
            kill_vector[h] = kill_tmp[acq->recpos_loc[3][i]][ishot];
            h++;
        }
    }
    /* ------------------------------- END TRACE KILLING ------------------------------- */

    /* Integration of measured and synthetic data  */
    for (i = 1; i <= gv->NTR; i++) {
        intseis_s = 0;
        intseis_sd = 0;
        if (vinv->VELOCITY == 0) {  /* only integration if displacement seismograms are inverted */
            for (j = 1; j <= gv->NS; j++) {
                intseis_s += section[i][j];
                intseis_section[i][j] = intseis_s * gv->DT;
                intseis_sd += sectiondata[i][j];
                intseis_sectiondata[i][j] = intseis_sd * gv->DT;
            }
        } else {
            for (j = 1; j <= gv->NS; j++) {
                intseis_section[i][j] = section[i][j];
                intseis_sectiondata[i][j] = sectiondata[i][j];
            }
        }
    }

    /* TIME WINDOWING */
    if (vinv->TIMEWIN == 1) {
        //time_window(intseis_section, iter, gv->NTRG,acq->recpos_loc, gv->NTR, gv->NS, ishot);
        //time_window(intseis_sectiondata, iter, gv->NTRG,acq->recpos_loc, gv->NTR, gv->NS, ishot);
    }

    /* NORMALIZE TRACES */
    if (vinv->NORMALIZE == 1) {
        //normalize_data(intseis_section,gv->NTR,gv->NS);
        //normalize_data(intseis_sectiondata,gv->NTR,gv->NS);
    }

    /* calculate RMS */
    /*for(i=1;i<=gv->NTR;i++){
     * for(j=1;j<=gv->NS;j++){ */
    /*RMS += (sectionpdata[i][j]-sectionvx[i][j]) * (sectionpdata[i][j]-sectionvx[i][j]); */
    /*RMS += (sectionpdata[i][j]) * (sectionpdata[i][j]); */
    /*      Lcount++;
     * }
     * } */

    /*RMS=sqrt(RMS/Lcount); */
    RMS = 1.0;

    /* calculate weighted data residuals and reverse time direction */
    /* calculate kind of "energy" */

    /* calculate envelope and hilbert transform for LNORM==8 */
    if (vinv->LNORM == 8) {
        //calc_envelope(intseis_sectiondata,intseis_sectiondata_envelope,gv->NS,gv->NTR);
        //calc_envelope(intseis_section,intseis_section_envelope,gv->NS,gv->NTR);
        //calc_hilbert(intseis_section,intseis_section_hilbert,gv->NS,gv->NTR);
        for (i = 1; i <= gv->NTR; i++) {
            for (j = 1; j <= gv->NS; j++) {
                dummy_1[i][j] =
                    (intseis_section_envelope[i][j] -
                     intseis_sectiondata_envelope[i][j]) / (intseis_section_envelope[i][j] +
                                                            vinv->WATERLEVEL_LNORM8) * intseis_section_hilbert[i][j];
            }
        }
        //calc_hilbert(dummy_1,dummy_2,gv->NS,gv->NTR);
    }
    /* end of calculate envelope and hilbert transform for LNORM==8 */

    for (i = 1; i <= gv->NTR; i++) {
        if ((vinv->TRKILL == 1) && (kill_vector[i] == 1))
            continue;

        if ((vinv->LNORM == 5) || (vinv->LNORM == 7)) {
            abs_sectiondata = 0.0;
            abs_section = 0.0;
            sectiondata_mult_section = 0.0;

            for (j = 1; j <= gv->NS; j++) {
                abs_sectiondata += intseis_sectiondata[i][j] * intseis_sectiondata[i][j];
                abs_section += intseis_section[i][j] * intseis_section[i][j];
                sectiondata_mult_section += intseis_sectiondata[i][j] * intseis_section[i][j];  /* calculation of dot product for measured (section) and synthetic (sectiondata) data */
            }
            if (abs_sectiondata == 0) {
                abs_sectiondata = 1;
            } else {
                abs_sectiondata = sqrt(abs_sectiondata);
            }
            if (abs_section == 0) {
                abs_section = 1;
            } else {
                abs_section = sqrt(abs_section);
            }
        }
        /* calculate residual seismograms and norm */

        /*if ((vinv->TRKILL == 1) && (kill_vector[i] == 1))
         * continue; */

        /*reverse time direction */
        invtime = gv->NS;

        for (j = 1; j <= gv->NS; j++) {
            /*printf("%d \t %d \t %e \t %e \n",i,j,sectionpdata[i][j],sectionp[i][j]); */
            /* calculate L1 residuals */
            if (vinv->LNORM == 1) {
                if (((sectiondata[i][j] - section[i][j]) / RMS) > 0)
                    signL1 = 1.0;
                if (((sectiondata[i][j] - section[i][j]) / RMS) < 0)
                    signL1 = -1.0;
                if (((sectiondata[i][j] - section[i][j]) / RMS) == 0)
                    signL1 = 0.0;

                sectiondiff[i][invtime] = signL1 / RMS;
            }

            /* calculate L2 residuals */
            if (vinv->LNORM == 2) {
                sectiondiff[i][invtime] = intseis_section[i][j] - intseis_sectiondata[i][j];
            }

            /* calculate Cauchy residuals *//* NOT UP TO DATE */
            if (vinv->LNORM == 3) {
                sectiondiff[i][invtime] =
                    ((sectiondata[i][j] - section[i][j]) / RMS) / (1 +
                                                                   (((sectiondata[i][j] -
                                                                      section[i][j]) / RMS) * ((sectiondata[i][j] -
                                                                                                section[i][j]) /
                                                                                               RMS))) / RMS;
            }

            /* calculate sech residuals *//* NOT UP TO DATE */
            if (vinv->LNORM == 4) {
                sectiondiff[i][invtime] = (tanh((sectiondata[i][j] - section[i][j]) / RMS)) / RMS;
            }

            /* calculate LNORM 5 or LNORM 7 residuals */
            if ((vinv->LNORM == 5) || (vinv->LNORM == 7)) {
                sectiondiff[i][invtime] =
                    ((intseis_section[i][j] * sectiondata_mult_section) /
                     (abs_section * abs_section * abs_section * abs_sectiondata)) -
                    (intseis_sectiondata[i][j] / (abs_section * abs_sectiondata));
            }

            /* calculate LNORM 8 residuals */
            if (vinv->LNORM == 8) {
                sectiondiff[i][invtime] =
                    intseis_section[i][j] / (intseis_section_envelope[i][j] +
                                             vinv->WATERLEVEL_LNORM8) * (intseis_section_envelope[i][j] -
                                                                         intseis_sectiondata_envelope[i][j]) -
                    dummy_2[i][j];
            }

            /*sectionpdiff[i][invtime]=sectionp[i][j]; */

            /* calculate norm for different LNORM */
            if ((vinv->LNORM == 2) && (swstestshot == 1)) {
                vinv->L2 += sectiondiff[i][invtime] * sectiondiff[i][invtime];
            }

            if (((vinv->LNORM == 5) || (vinv->LNORM == 7)) && (swstestshot == 1)) {
                vinv->L2 -= (intseis_sectiondata[i][j] * intseis_section[i][j]) / (abs_sectiondata * abs_section);
            }

            if ((vinv->LNORM == 8) && (swstestshot == 1)) {
                vinv->L2 +=
                    (intseis_section_envelope[i][invtime] -
                     intseis_sectiondata_envelope[i][invtime]) * (intseis_section_envelope[i][invtime] -
                                                                  intseis_sectiondata_envelope[i][invtime]);
            }

            if ((sws == 2) && (swstestshot == 1)) {
                vinv->L2 += fabs(sectiondiff[i][invtime]) * fabs(sectiondiffold[i][j]);
            }

            /*vinv->L2+=sectiondiff[i][invtime]; */

            invtime--;          /* reverse time direction */
        }
    }

    // taper(sectiondiff,gv->NTR,gv->NS);
    /* printf("\n MYID = %i   IN CALC_RES: L2 = %10.12f \n",MYID,l2); */

    /* free memory for integrated seismograms */
    free_matrix(intseis_section, 1, gv->NTR, 1, gv->NS);
    free_matrix(intseis_sectiondata, 1, gv->NTR, 1, gv->NS);

    /* free memory for time windowing and trace killing */
    if (vinv->TIMEWIN == 1)
        free_vector(picked_times, 1, gv->NTR);

    if (vinv->TRKILL == 1) {
        free_imatrix(kill_tmp, 1, gv->NTRG, 1, acq->nsrc);
        free_ivector(kill_vector, 1, gv->NTR);
    }

    if (vinv->LNORM == 8) {
        free_matrix(intseis_sectiondata_envelope, 1, gv->NTR, 1, gv->NS);
        free_matrix(intseis_section_envelope, 1, gv->NTR, 1, gv->NS);
        free_matrix(intseis_section_hilbert, 1, gv->NTR, 1, gv->NS);
        free_matrix(dummy_1, 1, gv->NTR, 1, gv->NS);
        free_matrix(dummy_2, 1, gv->NTR, 1, gv->NS);
    }
}
