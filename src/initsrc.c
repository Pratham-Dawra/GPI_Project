
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
 * Setting source parameters.
 * -------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "macros.h"
#include <stdio.h>

void initsrc(int ishot, int nshots, AcqVar *acq, GlobVar *gv)
{
    int nt;
    char sigf[STRING_SIZE * 2], file_ext[5];

    if (acq->nsrc_loc > 0) {
        free_matrix(acq->signals, 1, acq->nsrc_loc, 1, gv->NT);
        acq->signals = NULL;
        free_matrix(acq->srcpos_loc, 1, NSPAR, 1, 1);
        acq->srcpos_loc = NULL;
        acq->nsrc_loc = 0;
    }

    if (gv->RUN_MULTIPLE_SHOTS) {
        log_infoc(0, "Starting simulation for shot %d of %d.\n", ishot, nshots);
        float **srcpos_current = matrix(1, NSPAR, 1, 1);
        for (nt = 1; nt <= NSPAR; nt++) {
            srcpos_current[nt][1] = acq->srcpos[nt][ishot];
        }
        /* find this single source positions on subdomains  */
        acq->srcpos_loc = splitsrc(srcpos_current, &acq->nsrc_loc, 1, gv);
        free_matrix(srcpos_current, 1, NSPAR, 1, 1);
    } else {
        /* distribute source positions on subdomains */
        acq->srcpos_loc = splitsrc(acq->srcpos, &acq->nsrc_loc, acq->nsrc, gv);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* calculate wavelet for source point(s) */
    if (acq->nsrc_loc > 0) {
        wavelet(acq, gv);
    }

    /* write source wavelet to file in subdomain that contains source */
    if (acq->nsrc_loc > 0 && 1 == gv->SIGOUT) {
        switch (gv->SIGOUT_FORMAT) {
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
        int **dummy = imatrix(1, 3, 1, 1);
        dummy[1][1] = (int)(acq->srcpos_loc[1][1]) - 1 + gv->GGRID[1];
        dummy[2][1] = (int)(acq->srcpos_loc[2][1]) - 1 + gv->GGRID[3];
        dummy[3][1] = 0;
        sprintf(sigf, "%s.shot%d.%s", gv->SIGOUT_FILE, ishot, file_ext);
        log_info("Writing source wavelet to file %s.\n", sigf);
        FILE *fid = fopen(sigf, "w");
        if (!fid) {
            log_warn("Could not open %s for writing - no signature written.\n", sigf);
        } else {
            outseis_glob(fid, acq->signals, dummy, 1, acq->srcpos, gv->NT, gv->SIGOUT_FORMAT, ishot, 0, gv);
        }
        free_imatrix(dummy, 1, 3, 1, 1);
    }
}
