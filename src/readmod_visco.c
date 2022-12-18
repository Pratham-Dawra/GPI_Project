
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

/* ---------------------------------------------------------------------------
 * Read viscoelastic model properties (vp,vs,rho,qp,qs) from files
 * -------------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_su.h"
#include <unistd.h>
#include <stdbool.h>

void readmod_visco(float **rho, float **pi, float **u, float **taus, float **taup, float *eta, GlobVar *gv)
{
    float muv, piv, *pts, ts, tp, sumu, sumpi;
    int ii, jj;
    size_t ny;
    char filename[STRING_SIZE + 16];
    bool b_issu = false;

    const char *model[] = { "P-wave velocity model", "S-wave velocity model", "density model",
        "Quality factor Qp model", "Quality factor Qs model" };
    const char *suffix[] = { "vp", "vs", "rho", "qp", "qs" };
    const char *suffix_su[] = { "vp.su", "vs.su", "rho.su", "qp.su", "qs.su" };

    enum PARAMOD {
        P_VP = 0,
        P_VS,
        P_RHO,
        P_QP,
        P_QS,
        NPARA
    };

    FILE *fp[NPARA];

    /* vector for maxwell bodies */
    pts = vector(1, gv->L);
    for (int l = 1; l <= gv->L; l++) {
        pts[l] = 1.0 / (2.0 * PI * gv->FL[l]);
        eta[l] = gv->DT / pts[l];
    }

    float ws = 2.0 * PI / gv->TS;

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
        ny = (size_t)(gv->NYG);
    }

    float **para = (float **)malloc2d(NPARA, ny, sizeof(float));

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        if (b_issu) {
            su_read_trace(fp[P_VP], 0, (unsigned short)ny, false, NULL, &(para[P_VP][0]));
            su_read_trace(fp[P_VS], 0, (unsigned short)ny, false, NULL, &(para[P_VS][0]));
            su_read_trace(fp[P_RHO], 0, (unsigned short)ny, false, NULL, &(para[P_RHO][0]));
            su_read_trace(fp[P_QP], 0, (unsigned short)ny, false, NULL, &(para[P_QP][0]));
            su_read_trace(fp[P_QS], 0, (unsigned short)ny, false, NULL, &(para[P_QS][0]));
        } else {
            fread(&(para[P_VP][0]), sizeof(float), ny, fp[P_VP]);
            fread(&(para[P_VS][0]), sizeof(float), ny, fp[P_VS]);
            fread(&(para[P_RHO][0]), sizeof(float), ny, fp[P_RHO]);
            fread(&(para[P_QP][0]), sizeof(float), ny, fp[P_QP]);
            fread(&(para[P_QS][0]), sizeof(float), ny, fp[P_QS]);
        }
        for (int j = 1; j <= gv->NYG; j++) {
            tp = 2.0 / (para[P_QP][j - 1] * gv->L);
            ts = 2.0 / (para[P_QS][j - 1] * gv->L);
            muv = para[P_VS][j - 1] * para[P_VS][j - 1] * para[P_RHO][j - 1];
            piv = para[P_VP][j - 1] * para[P_VP][j - 1] * para[P_RHO][j - 1];
            sumu = 0.0f;
            sumpi = 0.0f;
            for (int l = 1; l <= gv->L; l++) {
                sumu += ((ws * ws * pts[l] * pts[l] * ts) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumpi += ((ws * ws * pts[l] * pts[l] * tp) / (1.0 + ws * ws * pts[l] * pts[l]));
            }
            muv = muv / (1.0 + sumu);
            piv = piv / (1.0 + sumpi);
            /* only the PE which belongs to the current global gridpoint 
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;
                taus[jj][ii] = ts;
                taup[jj][ii] = tp;
                u[jj][ii] = muv;
                rho[jj][ii] = para[P_RHO][j - 1];
                pi[jj][ii] = piv;
            }
        }
    }

    free_vector(pts, 1, gv->L);
    if (para)
        free(para);

    for (int i = 0; i < NPARA; ++i) {
        fclose(fp[i]);
    }

    return;
}
