
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

int main(int argc, char **argv)
{
    int nsnap;
    GlobVar gv = {.MPID = 0,.OUTNTIMESTEPINFO = 1,.NDT = 1,.IDX = 1,.IDY = 1 };

    log_init(NULL);
    log_banner(LOG_SNAP);

    if (argc != 2) {
        log_fatal
            ("Unexpected number of commandline arguments; single argument required: name of json parameter file.\n");
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

    gv.NXG = gv.NX;
    gv.NYG = gv.NY;
    gv.NX = gv.NXG / gv.NPROCX;
    gv.NY = gv.NYG / gv.NPROCY;

    nsnap = 1 + iround((gv.TSNAP2 - gv.TSNAP1) / gv.TSNAPINC);

    switch (gv.SNAP) {
      case 1:                  /*particle velocity */
          merge(nsnap, 1, &gv);
          merge(nsnap, 2, &gv);
          break;
      case 2:                  /*pressure */
          merge(nsnap, 6, &gv);
          break;
      case 4:                  /*particle velocity */
          merge(nsnap, 1, &gv);
          merge(nsnap, 2, &gv);
          merge(nsnap, 6, &gv);
          /* FALL THRU */
      case 3:
          merge(nsnap, 4, &gv);
          merge(nsnap, 5, &gv);
          break;
      default:
          log_fatal("Unknown value for parameter SNAP.\n");
          break;
    }

    log_finalize();

    return EXIT_SUCCESS;
}
