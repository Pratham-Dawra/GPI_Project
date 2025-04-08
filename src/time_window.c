
/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * --------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Apply time damping (after Brossier (2009))                                 
 *  ----------------------------------------------------------------------*/

#include <stdio.h>

#include "fd.h"
#include "logging.h"
#include "util.h"

void time_window(float **sectiondata1, float **sectiondata2, int **recpos_loc, int ntr,
                 int ntr_loc, int ishot, int groupnum, int stf_switch, const GlobVar *gv){

    /* declaration of local variables */
    char pickfile[2*STRING_SIZE], dummy1[STRING_SIZE], dummy2[STRING_SIZE];
    float time, dumpa, dumpb;
    float taper_top = 0.0, taper_base = 0.0;
    int i, j;

    FILE *fptime = NULL;

    if (ntr_loc == 0)
        return;

    float **picked_times_m = fmatrix(1,2,1,ntr_loc);
    float **pick_tmp_m = fmatrix(1,2,1,ntr);

    /* read picked first arrival times */
    if(stf_switch){
        sprintf(pickfile,"%s_shot%d_stage%d_stf.dat",gv->TW_PICKFILE,ishot,groupnum);
        fptime=fopen(pickfile,"r");
        if (fptime == NULL) {
            sprintf(pickfile,"%s_shot%d_stf.dat",gv->TW_PICKFILE,ishot);
            fptime=fopen(pickfile,"r");
            if (fptime == NULL) {
                sprintf(pickfile,"%s_shot%d.dat",gv->TW_PICKFILE,ishot);
                fptime=fopen(pickfile,"r");
                if (fptime == NULL) {
                    log_fatal("%s_shot%d.dat could not be opened!\n", gv->TW_PICKFILE,ishot);
                }
            }
        }
    }else{
        sprintf(pickfile,"%s_shot%d_stage%d.dat",gv->TW_PICKFILE,ishot,groupnum);
        fptime=fopen(pickfile,"r");
        if (fptime == NULL) {
            sprintf(pickfile,"%s_shot%d.dat",gv->TW_PICKFILE,ishot);
            fptime=fopen(pickfile,"r");
            if (fptime == NULL) {
                log_fatal("%s_shot%d.dat could not be opened!\n", gv->TW_PICKFILE,ishot);
            }
        }
    }

    //log_info("Read mute from pick file %s, starting at line %d.\n", pickfile, gv->TW_START_READING_AT_LINE);

    for(i=1;i<=gv->TW_START_READING_AT_LINE-1;i++){
        fscanf(fptime,"%s %s\n", dummy1, dummy2);
    }
    for(i=1;i<=ntr;i++){
        fscanf(fptime,"%f %f",&pick_tmp_m[1][i],&pick_tmp_m[2][i]);
    }

    fclose(fptime);

    /* distribute picks on CPUs */
    for(i=1;i<=ntr_loc;i++){
        picked_times_m[1][i] = pick_tmp_m[1][recpos_loc[4][i]];
        picked_times_m[2][i] = pick_tmp_m[2][recpos_loc[4][i]];
    }

    free_matrix(pick_tmp_m,1,2,1,ntr);

    for(i=1;i<=ntr_loc;i++){
        for(j=2;j<=gv->NT;j++){
            time = (float)(j * gv->DT);
            
            if (1 == gv->USE_TW) {
                dumpa = (time-picked_times_m[1][i]);
                taper_top = exp(gv->TW_EXP_START*dumpa);

                dumpb = (time-picked_times_m[2][i]); 
                taper_base = exp(-gv->TW_EXP_END*dumpb);

                if(time<=picked_times_m[1][i]){
                    sectiondata1[i][j] = sectiondata1[i][j] * taper_top;
                    sectiondata2[i][j] = sectiondata2[i][j] * taper_top;
                }
                if(time>=picked_times_m[2][i]){
                    sectiondata1[i][j] = sectiondata1[i][j] * taper_base;
                    sectiondata2[i][j] = sectiondata2[i][j] * taper_base;
                }
            } else if (2 == gv->USE_TW) {
                if (gv->TW_SIN_START >= 0) {
                    if (time < (picked_times_m[1][i] - gv->TW_SIN_START)) {
                        taper_top = 0;
                    } else if (time > picked_times_m[1][i]) {
                        taper_top = 1;
                    } else {
                        dumpa = (picked_times_m[1][i] - time) / gv->TW_SIN_START;
                        taper_top = 0.5 * (1 + cosf(dumpa * M_PI));
                    }
                } else {
                    if (time < picked_times_m[1][i]) {
                        taper_top = 0;
                    } else if (time > (picked_times_m[1][i] - gv->TW_SIN_START)) {
                        taper_top = 1;
                    } else {
                        dumpa = (picked_times_m[1][i] - gv->TW_SIN_START - time) / gv->TW_SIN_START;
                        taper_top = 0.5 * (1 + cosf(dumpa * M_PI));
                    }
                }
                if (gv->TW_SIN_END >= 0) {
                    if (time < (picked_times_m[2][i])) {
                        taper_base = 1;
                    } else if (time > (picked_times_m[2][i] + gv->TW_SIN_END)) {
                        taper_base = 0;
                    } else {
                        dumpb = (picked_times_m[2][i] - time) / gv->TW_SIN_END;
                        taper_base = 0.5 * (1 + cosf(dumpb * M_PI));
                    }
                } else {
                    if (time < (picked_times_m[2][i] + gv->TW_SIN_END)) {
                        taper_base = 1;
                    } else if (time > picked_times_m[2][i]) {
                        taper_base = 0;
                    } else {
                        dumpb = (picked_times_m[2][i] - gv->TW_SIN_END - time) / gv->TW_SIN_END;
                        taper_base = 0.5 * (1 + cosf(dumpb * M_PI));
                    }
                }
                sectiondata1[i][j] = sectiondata1[i][j] * taper_top;
                sectiondata2[i][j] = sectiondata2[i][j] * taper_top;
                sectiondata1[i][j] = sectiondata1[i][j] * taper_base;
                sectiondata2[i][j] = sectiondata2[i][j] * taper_base;
                
            }
        }
    }

    free_matrix(picked_times_m,1,2,1,ntr_loc);

} /* end of function time_window.c */
