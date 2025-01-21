
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

/* ------------------------------------------------------------------------
 * Scan parameter file for topography to later place sources and receivers
 * at a specific point below the topography.
 * -----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_su.h"
#include <unistd.h>
#include <stdbool.h>
#include <float.h>

int *scan_topo(GlobVar *gv)
{
    const char *model = "P-wave velocity model";
    const char *suffix = "vp";
    const char *suffix_su = "vp.su";

    char filename[STRING_SIZE + 16];
    bool b_issu = false;

    /* test for SU model - otherwise we read standard binary */
    sprintf(filename, "%s.%s", gv->MFILE, suffix_su);
    if (!access(filename, R_OK)) {
        b_issu = true;
    } else {
        sprintf(filename, "%s.%s", gv->MFILE, suffix);
        if (access(filename, R_OK)) {
            log_fatal("Could not open %s to scan topography.\n", model);
        }
    }

    FILE *fp = fopen(filename, "rb");
    if (!fp)
        log_fatal("Could not open %s %s.\n", model, filename);

    int *topo = ivector(1, gv->NXG);
    float *vptrace = (float *)malloc1d(gv->NYG, sizeof(float));
    float vp_prev;

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        if (b_issu) {
            su_read_trace(fp, 0, (unsigned short)gv->NYG, false, NULL, &(vptrace[0]));
        } else {
            fread(&(vptrace[0]), sizeof(float), gv->NYG, fp);
        }
        topo[i] = 1;
        vp_prev = vptrace[0];
        for (int j = 2; j <= gv->NYG; j++) {
            if (vptrace[j - 1] != vp_prev) {
                topo[i] = j;
                break;
            }
        }
        log_debug("Topography scan: trace %d, topo grid point %d (depth: %f)\n", i, topo[i], (topo[i] - 1) * gv->DH);
    }

    free(vptrace);
    fclose(fp);

    return topo;
}
