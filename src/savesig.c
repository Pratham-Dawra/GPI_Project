
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
--------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * write seismograms to files
 *----------------------------------------------------------------------*/

#include <stdio.h>

#include "fd.h"
#include "util.h"
#include "logging.h"
#include "macros.h"
#include "su_struct.h"
#include "read_su.h"
#include "write_su.h"

void savesig(float **signals, AcqVar *acq, int nsrc_loc, int ishot, int iter, int fswitch, GlobVar *gv)
{
    int nsig = 0;
    char sigf[2 * STRING_SIZE], seisfile[STRING_SIZE], file_ext[8], file_iter[16] = "";
    FILE *fp = NULL;

    float **srcpos1 = fmatrix(1, 7, 1, 1);
    for (int nt = 1; nt <= 7; nt++) {
        srcpos1[nt][1] = acq->srcpos[nt][ishot];
    }

    switch (gv->SIGOUT_FORMAT) {
      case 0:
          /* FALLTHRU */
      case 1:
          sprintf(file_ext, "su");
          break;
      case 2:
          sprintf(file_ext, "txt");
          break;
      case 3:
          sprintf(file_ext, "bin");
          break;
    }

    sprintf(seisfile, "%s", gv->SIGOUT_FILE);

    if (gv->METHOD) {
        sprintf(file_iter, "%s%d", "_iter", iter);
    }

    switch (fswitch) {
      case 0:
          if (gv->RUN_MULTIPLE_SHOTS) {
              sprintf(sigf, "%s_srcsig%s.%s.shot%d", seisfile, file_iter, file_ext, ishot);
              nsig = 1;
          } else {
              sprintf(sigf, "%s_srcsig%s.%s", seisfile, file_iter, file_ext);
              nsig = acq->nsrc;
          }
          break;
      case 1:
          if (gv->RUN_MULTIPLE_SHOTS) {
              sprintf(sigf, "%s_srcsig_stfi%s.%s.shot%d", seisfile, file_iter, file_ext, ishot);
              nsig = 1;
          } else {
              sprintf(sigf, "%s_srcsig_stfi%s.%s", seisfile, file_iter, file_ext);
              nsig = acq->nsrc;
          }
          break;
      case 2:
          if (gv->RUN_MULTIPLE_SHOTS) {
              sprintf(sigf, "%s_srcsig_stfi_stepl%s.%s.shot%d", seisfile, file_iter, file_ext, ishot);
              nsig = 1;
          } else {
              sprintf(sigf, "%s_srcsig_stfi_stepl%s.%s", seisfile, file_iter, file_ext);
              nsig = acq->nsrc;
          }
          break;
    }

    if (gv->RUN_MULTIPLE_SHOTS) {
        log_warn("value of nsrc loc is %d in if \n", nsrc_loc);
        if (nsrc_loc) {
            fp = fopen(sigf, "w");
            if (!fp) {
                log_error("Could not open %s for writing. Continuing.\n", sigf);
            } else {
                outseis(fp, gv->SOURCE_TYPE, signals, acq->recpos, NULL, nsig, srcpos1, 1, ishot, gv);
                log_info("Writing source signature to file %s.\n", sigf);
            }
        }
    } else {
         log_warn("value of nsrc loc is %d in else \n", nsrc_loc);
        float **fulldata = matrix(1, nsig, 0, gv->NS);
        log_info("fulldata is %");
        fp = fopen(sigf, "w");
        if (0 == gv->MPID) {
            if (!fp) {
                log_error("Could not open %s for writing. Continuing.\n", sigf);
        } else {
            catseis(signals, fulldata, acq->srcswitch, nsig, gv->NS);
            outseis(fp, gv->SOURCE_TYPE, fulldata, acq->recpos, NULL, nsig, srcpos1, 1, ishot, gv);
            log_info("Writing source signature to file %s.\n", sigf);
            }
        }

            free_matrix(fulldata, 1, nsig, 0, gv->NS-1);

    }
    free_matrix(srcpos1, 1, 7, 1, 1);
    return;
}
