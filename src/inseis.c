
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
 *   Write seismograms to disk                                  
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "logging.h"
//#include "segy.h"
#include "read_su.h"

void inseis(int ishot, int sws, float **section, GlobVar *gv, GlobVarInv *vinv)
{
    /* declaration of extern variables */
    char data[STRING_SIZE * 2];
    FILE *fpdata;
    int i;

    if (sws == 1) {             /* open seismic data vx */
        sprintf(data, "%s_vx.su.shot%d", vinv->DATA_DIR, ishot);
        if (gv->MPID == 0)
            log_info("Read: %s\n", data);
    }

    if (sws == 2) {             /* open seismic data vy */
        sprintf(data, "%s_vy.su.shot%d", vinv->DATA_DIR, ishot);
        if (gv->MPID == 0)
            log_info("Read: %s\n", data);
    }

    /* sws 3 -- 8 not used */

    //if(sws==7){  /* open convolved seismic data vx */
    //   sprintf(data,"%s_vx.su.conv.shot%d",vinv->DATA_DIR,ishot);
    //}

    //if(sws==8){  /* open convolved seismic data vy */
    //    sprintf(data,"%s_vy.su.conv.shot%d",vinv->DATA_DIR,ishot);
    //}

    if (sws == 9) {             /* open seismic data p */
        sprintf(data, "%s_p.su.shot%d", vinv->DATA_DIR, ishot);
        if (gv->MPID == 0)
            log_info("Read: %s\n", data);
    }
    //if(sws==10){  /* open seismic data vz */
    //    sprintf(data,"%s_vz.su.shot%d",vinv->DATA_DIR,ishot);
    //}

    fpdata = fopen(data, "r");
    if (fpdata == NULL) {
        if (gv->MPID == 0)
            log_info("Was not able to read %s\n", data);
        log_fatal("Seismograms for inversion were not found.\n");
    }

    /* Read su data */
    for (i = 1; i <= gv->NTRG; i++) {
        su_read_trace(fpdata, 0, (unsigned short)gv->NS, false, NULL, &(section[i][1]));
    }
    fclose(fpdata);
}
