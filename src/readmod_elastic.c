
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
 * Read elastic model properties (vp,vs,rho) from files
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_su.h"
#include <unistd.h>
#include <stdbool.h>
#include <float.h>

void readmod_elastic(MemModel *mpm, MemInv *minv, GlobVar *gv)
{
    int ii, jj;
    size_t ny;
    char filename[STRING_SIZE + 16];
    bool b_issu = false;

    const char *model[] = { "P-wave velocity model", "S-wave velocity model", "density model" };
    const char *suffix[] = { "vp", "vs", "rho" };
    const char *suffix_su[] = { "vp.su", "vs.su", "rho.su" };

    enum PARAMOD {
        P_VP = 0,
        P_VS,
        P_RHO,
        NPARA
    };

    FILE *fp[NPARA];

    /* test for SU model - otherwise we read standard binary */
    sprintf(filename, "%s.%s", gv->MFILE, suffix_su[0]);
    if (!access(filename, R_OK)) {
        b_issu = true;
        log_infoc(0, "Reading models in SU format.\n");
        sprintf(filename, "%s.%s", gv->MFILE, suffix[0]);
        if (!access(filename, R_OK)) {
            log_warnc(0, "Models are also available as plain binary files. Check consistency, if applicable.\n");
        }
    }

    for (int i = 0; i < NPARA; ++i) {
        if (b_issu) {
            sprintf(filename, "%s.%s", gv->MFILE, suffix_su[i]);
        } else {
            sprintf(filename, "%s.%s", gv->MFILE, suffix[i]);
        }
        log_infoc(0, "Reading %s: %s\n", model[i], filename);
        fp[i] = fopen(filename, "rb");
        if (!fp[i])
            log_fatal("Could not open %s %s.\n", model[i], filename);
    }

    /* if we read SU format, we can perform additional cross-checks; however, coordinates are not used atm */
    if (b_issu) {
        size_t nx;
        unsigned short ns, dt;
        short int delrt;
        for (int i = 0; i < NPARA; ++i) {
            nx = su_get_nt(fp[i], &ns, &dt, &delrt);
            if (nx < (unsigned)gv->NXG)
                log_fatal("%s has fewer than NX=%d traces.\n", model[i], gv->NXG);
            else if (nx > (unsigned)gv->NXG)
                log_warnc(0, "%s has more than NX=%d traces; ignoring add. traces.\n", model[i], gv->NXG);
            if (ns < (unsigned short)gv->NYG)
                log_fatal("%s has fewer than NY=%d samples.\n", model[i], gv->NYG);
            else if (ns > (unsigned short)gv->NYG)
                log_warnc(0, "%s has more than NY=%d samples; ignoring add. samples.\n", model[i], gv->NYG);
        }
        ny = ns;
    } else {
        ny = (size_t) (gv->NYG);
    }

    float **para = (float **)malloc2d(NPARA, ny, sizeof(float));
    gv->VPMIN = FLT_MAX;
    gv->VPMAX = 0.0;
    gv->VSMIN = FLT_MAX;
    gv->VSMAX = 0.0;

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        if (b_issu) {
            su_read_trace(fp[P_VP], 0, (unsigned short)ny, false, NULL, &(para[P_VP][0]));
            su_read_trace(fp[P_VS], 0, (unsigned short)ny, false, NULL, &(para[P_VS][0]));
            su_read_trace(fp[P_RHO], 0, (unsigned short)ny, false, NULL, &(para[P_RHO][0]));
        } else {
            fread(&(para[P_VP][0]), sizeof(float), ny, fp[P_VP]);
            fread(&(para[P_VS][0]), sizeof(float), ny, fp[P_VS]);
            fread(&(para[P_RHO][0]), sizeof(float), ny, fp[P_RHO]);
        }
        for (int j = 1; j <= gv->NYG; j++) {
            float vp = para[P_VP][j-1];
            float vs = para[P_VS][j-1];
	        if (vp < gv->VPMIN && vp > V_IGNORE) gv->VPMIN = vp;
	        if (vp > gv->VPMAX && vp > V_IGNORE) gv->VPMAX = vp;
	        if (vs < gv->VSMIN && vs > V_IGNORE) gv->VSMIN = vs;
	        if (vs > gv->VSMAX && vs > V_IGNORE) gv->VSMAX = vs;
            /* only the PE which belongs to the current global gridpoint 
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;
                if (gv->MODE == FWI) {
                    minv->Vp0[jj][ii] = vp;
                    minv->Vs0[jj][ii] = vs;
                    minv->Rho0[jj][ii] = para[P_RHO][j - 1];
                }
                mpm->pu[jj][ii] = vs * vs * para[P_RHO][j - 1];
                mpm->prho[jj][ii] = para[P_RHO][j - 1];
                mpm->ppi[jj][ii] = vp * vp * para[P_RHO][j - 1];
            }
        }
    }

    if (para)
        free(para);

    for (int i = 0; i < NPARA; ++i) {
        fclose(fp[i]);
    }

    return;
}
