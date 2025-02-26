
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

/*----------------------------------------------------------------------
 * merge snapshots files written by the different processes to
 * a single file
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include "fd.h"
#include "logging.h"

void merge(int nsnap, int type, int SNAPIDX[][5], GlobVar *gv)
{
    char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10];
    FILE *fp[gv->NPROCY][gv->NPROCX], *fpout = NULL;
    float a;

    switch ((gv->SNAP_FORMAT)) {
      case 1:
          sprintf(ext, ".su");
          break;
      case 2:
          sprintf(ext, ".asc");
          break;
      case 3:
          sprintf(ext, ".bin");
          break;
    }

    switch (type) {
      case 1:
          log_info("Merging x-component of particle velocity");
          strcat(ext, ".vx");
          break;
      case 2:
          log_info("Merging y-component of particle velocity");
          strcat(ext, ".vy");
          break;
      case 4:
          log_info("Merging P-wave energyfield");
          strcat(ext, ".div");
          break;
      case 5:
          log_info("Merging S-wave energyfield");
          strcat(ext, ".curl");
          break;
      case 6:
          log_info("Merging pressure");
          strcat(ext, ".p");
          break;
      default:
          log_fatal("Unknown type of snapshot files to merge.\n");
          break;
    }

    sprintf(mfile, "%s%s", (gv->SNAP_FILE), ext);
    log_std(" (files: %s.*.*).\n", mfile);

    sprintf(outfile, "%s%s", (gv->SNAP_FILE), ext);

    log_info("Output file: %s\n", outfile);

    fpout = fopen(outfile, "w");
    log_debug("Opening snapshot files %s.*.* for reading.\n", mfile);

    for (int ip = 0; ip < gv->NPROCX; ip++) {
        for (int jp = 0; jp < gv->NPROCY; jp++) {
            sprintf(file, "%s.%i.%i", mfile, ip, jp);
            fp[jp][ip] = fopen(file, "r");
            if (fp[jp][ip] == NULL)
                log_fatal("Could not open file %s for reading.\n", file);
        }
    }
    log_debug("Finished opening all input files. Now copying the data.\n");

    float max = 0.0;
    int n1 = 0, n2 = 0;
    for (int n = 0; n < nsnap; ++n) {
        for (int ip = 0; ip < gv->NPROCX; ip++) {
            for (int i = SNAPIDX[ip][1]; i <= SNAPIDX[ip][2]; i += gv->IDX) {
                if (0 == n) ++n2;
                for (int jp = 0; jp < gv->NPROCY; jp++) {
                    for (int j = SNAPIDX[jp*gv->NPROCX][3]; j <= SNAPIDX[jp*gv->NPROCX][4]; j += gv->IDY) {
                        a = readdsk(fp[jp][ip], (gv->SNAP_FORMAT));
                        if (a > max) max = a;
                        writedsk(fpout, a, gv->SNAP_FORMAT);
                        if (0 == n && 1 == n2) ++n1;
                    }
                }
            }
        }
    }

    log_debug("Finished copying the data. Now closing all input files.\n");

    for (int ip = 0; ip < gv->NPROCX; ip++) {
        for (int jp = 0; jp < gv->NPROCY; jp++) {
            fclose(fp[jp][ip]);
        }
    }
    log_debug("Finished closing all input files.\n");

    if ((gv->SNAP_FORMAT) == 3) {
        log_info("Use...\n");
        log_std
            ("%sxmovie n1=%d n2=%d < %s loop=1 label1=\"Y\" label2=\"X\" title=\"%%g\" d1=%f d2=%f f1=%f f2=%f clip=%e sleep=1%s\n",
             LOG_COLOR_BOLD, n1, n2, outfile, gv->DH*gv->IDY, gv->DH*gv->IDX, 0.0, 0.0, max / 10.0, LOG_ALL_RESET);
        log_info("...to play the snapshot movie.\n");
    }
}
