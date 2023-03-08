
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

#include <unistd.h>
#include "fd.h"
#include "logging.h"
#include "acq_struct.h"

int main(int argc, char **argv)
{
    GlobVar gv = {.MPID = 0,.OUTNTIMESTEPINFO = 1,.NDT = 1,.IDX = 1,.IDY = 1 };
    AcqVar acq = { };

    log_init(NULL);
    log_banner(LOG_SNAP);

    if (argc != 2) {
        log_fatal("Unexpected number of commandline arguments; single argument required: name of json parameter file.\n");
    }

    const char *fileinp = argv[1];

    log_info("Parameter json file: %s\n", fileinp);

    if (access(fileinp, R_OK)) {
        log_fatal("Could not open file %s for reading.\n", fileinp);
    }

    if (strstr(fileinp, ".json")) {
        read_par_json(fileinp, &gv);
    } else {
        log_fatal("Parameter file has no json suffix.\n");
    }

    // number of subdomains
    gv.NPROC = gv.NPROCX * gv.NPROCY;

    // domain size
    gv.NX = gv.NXG / gv.NPROCX;
    gv.NY = gv.NYG / gv.NPROCY;

    // number of snapshots (simulate same "algorithm" as used in SOFI)
    int lsnap = iround(gv.TSNAP1 / gv.DT);      /* first snapshot at this time step */
    int isnap = iround(gv.TSNAPINC / gv.DT);    /* snapshot increment in number of timesteps */
    int esnap = iround(gv.TSNAP2 / gv.DT);      /* last snapshot no later than this time step */
    int nsnap = 0;
    for (int nt = 1; nt <= gv.NT; nt++) {
        if ((nt == lsnap) && (nt <= esnap)) {
            ++nsnap;
            lsnap += isnap;
        }
    }
    log_info("Number of snapshots per shot: %d\n", nsnap);

    int nshots = 1;
    sources(&acq, &gv);
    if (gv.RUN_MULTIPLE_SHOTS) {
      nshots = acq.nsrc;
      log_info("Number of individual shots: %d\n", nshots);
    } else {
      log_info("All shots were run simultaneously.\n");
    }
    nsnap *= nshots;
    log_info("Total number of snapshots: %d\n", nsnap);

    // simulate MPI setup of SOFI without actually using MPI; see initproc.c for details
    int POS[gv.NPROC][3];
    int GGRID[gv.NPROC][5];
    int SNAPIDX[gv.NPROC][5];

    for (int irank = 0; irank < gv.NPROC; ++irank) {
        POS[irank][1] = irank % gv.NPROCX;
        POS[irank][2] = irank / gv.NPROCX;
        GGRID[irank][1] = 1;
        for (int i = 0; i < POS[irank][1]; ++i) {
            GGRID[irank][1] += gv.NX;
        }
        GGRID[irank][2] = GGRID[irank][1] + (gv.NX - 1);
        GGRID[irank][3] = 1;
        for (int i = 0; i < POS[irank][2]; ++i) {
            GGRID[irank][3] += gv.NY;
        }
        GGRID[irank][4] = GGRID[irank][3] + (gv.NY - 1);

        log_debug("PE %d; col %d; row %d; grid points xb-xe %d,%d yb-ye %d,%d\n",
                  irank, POS[irank][1], POS[irank][2],
                  GGRID[irank][1], GGRID[irank][2], GGRID[irank][3], GGRID[irank][4]);

        if (gv.IDX == 1) {
            SNAPIDX[irank][1] = GGRID[irank][1];
            SNAPIDX[irank][2] = GGRID[irank][2];
        } else {
            SNAPIDX[irank][1] = gv.NXG+1;
            SNAPIDX[irank][2] = -1;
            for (int i = 1; i <= gv.NXG; i += gv.IDX) {
                if (i >= GGRID[irank][1] && i <= GGRID[irank][2]) {
                    if (i < SNAPIDX[irank][1]) SNAPIDX[irank][1] = i;
                    if (i > SNAPIDX[irank][2]) SNAPIDX[irank][2] = i;
                }
            }
        }
        if (gv.IDY == 1) {
            SNAPIDX[irank][3] = GGRID[irank][3];
            SNAPIDX[irank][4] = GGRID[irank][4];
        } else {
            SNAPIDX[irank][3] = gv.NYG+1;
            SNAPIDX[irank][4] = -1;
            for (int j = 1; j <= gv.NYG; j += gv.IDY) {
                if (j >= GGRID[irank][3] && j <= GGRID[irank][4]) {
                    if (j < SNAPIDX[irank][3]) SNAPIDX[irank][3] = j;
                    if (j > SNAPIDX[irank][4]) SNAPIDX[irank][4] = j;
                }
            }
        }
        SNAPIDX[irank][1] = SNAPIDX[irank][1] - GGRID[irank][1] + 1;
        SNAPIDX[irank][2] = SNAPIDX[irank][2] - GGRID[irank][1] + 1;
        SNAPIDX[irank][3] = SNAPIDX[irank][3] - GGRID[irank][3] + 1;
        SNAPIDX[irank][4] = SNAPIDX[irank][4] - GGRID[irank][3] + 1;

        log_debug("PE %d; col %d; row %d; snap indices xs-xe %d,%d yb-ye %d,%d\n",
                  irank, POS[irank][1], POS[irank][2], 
                  SNAPIDX[irank][1], SNAPIDX[irank][2], SNAPIDX[irank][3], SNAPIDX[irank][4]);
    }

    switch (gv.SNAP) {
      case 1:                  /*particle velocity */
          merge(nsnap, 1, SNAPIDX, &gv);
          merge(nsnap, 2, SNAPIDX, &gv);
          break;
      case 2:                  /*pressure */
          merge(nsnap, 6, SNAPIDX, &gv);
          break;
      case 4:                  /*particle velocity */
          merge(nsnap, 1, SNAPIDX, &gv);
          merge(nsnap, 2, SNAPIDX, &gv);
          merge(nsnap, 6, SNAPIDX, &gv);
          /* FALL THRU */
      case 3:
          merge(nsnap, 4, SNAPIDX, &gv);
          merge(nsnap, 5, SNAPIDX, &gv);
          break;
      default:
          log_fatal("Unknown value for parameter SNAP.\n");
          break;
    }

    log_finalize();

    return EXIT_SUCCESS;
}
