
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
 * Read elastic TTI model properties (vp,vs,rho,eps,delta,theta) from files
 * -------------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_su.h"
#include <unistd.h>
#include <stdbool.h>
#include <float.h>

void readmod_acoustic_tti(MemModel * mpm, GlobVar * gv)
{
    float c11, c33, c13, c55, t;
    int ii, jj;
    size_t ny;
    float l1, l2, l12, l22, l14, l24;
    float a1, a3;
    char filename[STRING_SIZE + 16];
    bool b_issu = false;

    const char *model[] = { "P-wave velocity model", "density model",
        "Epsilon model", "Delta model", "Theta model"
    };
    const char *suffix[] = { "vp", "rho", "epsilon", "delta", "theta" };
    const char *suffix_su[] = { "vp.su", "rho.su", "epsilon.su", "delta.su", "theta.su" };

    enum PARAMOD {
        P_VP = 0,
        P_RHO,
        P_EPS,
        P_DEL,
        P_TET,
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
    gv->VSMIN = 0.0;
    gv->VSMAX = 0.0;

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        if (b_issu) {
            su_read_trace(fp[P_VP], 0, (unsigned short)ny, false, NULL, &(para[P_VP][0]));
            su_read_trace(fp[P_RHO], 0, (unsigned short)ny, false, NULL, &(para[P_RHO][0]));
            su_read_trace(fp[P_EPS], 0, (unsigned short)ny, false, NULL, &(para[P_EPS][0]));
            su_read_trace(fp[P_DEL], 0, (unsigned short)ny, false, NULL, &(para[P_DEL][0]));
            su_read_trace(fp[P_TET], 0, (unsigned short)ny, false, NULL, &(para[P_TET][0]));
        } else {
            fread(&(para[P_VP][0]), sizeof(float), ny, fp[P_VP]);
            fread(&(para[P_RHO][0]), sizeof(float), ny, fp[P_RHO]);
            fread(&(para[P_EPS][0]), sizeof(float), ny, fp[P_EPS]);
            fread(&(para[P_DEL][0]), sizeof(float), ny, fp[P_DEL]);
            fread(&(para[P_TET][0]), sizeof(float), ny, fp[P_TET]);
        }
        for (int j = 1; j <= gv->NYG; j++) {
            float vp = para[P_VP][j-1];
            float vp_90 = vp * sqrt(2.0 * para[P_EPS][j - 1] + 1.0);
            // epsilon can (theoretically) be negative, delta as well; make sure we get correct min/max
	        float vpmax = vp_90 > vp ? vp_90 : vp;
            float vpmin = vp_90 < vp ? vp_90 : vp;
            if (vpmin < gv->VPMIN && vp > V_IGNORE) gv->VPMIN = vpmin;
            if (vpmax > gv->VPMAX && vp > V_IGNORE) gv->VPMAX = vpmax;
            /* only the PE which belongs to the current global gridpoint 
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;
                c33 = para[P_RHO][j - 1] * vp * vp;
                c11 = c33 * (2.0 * para[P_EPS][j - 1] + 1.0);
                c55=0.0;
                c13 = sqrt((2.0 * para[P_DEL][j - 1] * c33 * (c33 - c55)) + ((c33 - c55) * (c33 - c55))) - c55;
                /* Bond transformation (Oh et al, 2020, GJI, doi: 10.1093/gji/ggaa295 */
                t = para[P_TET][j - 1] * PI / 180.0;
                l1 = cos(t);
                l2 = sin(t);
                l12 = l1 * l1;
                l22 = l2 * l2;
                l14 = l12 * l12;
                l24 = l22 * l22;
                /*l13 = l1 * l12;
                l23 = l2 * l22;*/
                a1 = 2.0 * c13 + 4.0 * c55;
                a3 = c11 + c33 - 4.0 * c55;
                /*a4 = c11 + c33 - 2.0 * c13;
                a5 = c13 - c11 + 2.0 * c55;
                a6 = c13 - c33 + 2.0 * c55;*/
                mpm->pc11[jj][ii] = c11 * l14 + c33 * l24 + a1 * l12 * l22;           // c11t
                mpm->prho[jj][ii] = para[P_RHO][j - 1];
                mpm->pc33[jj][ii] = c11 * l24 + c33 * l14 + a1 * l12 * l22;           // c33t
                mpm->pc13[jj][ii] = a3 * l12 * l22 + c13 * (l14 + l24);               // c13t
            }
        }
    }

    if (para)
        free(para);

    for (int i = 0; i < NPARA; ++i) {
        fclose(fp[i]);
    }

#ifdef EBUG
    debug_check_matrix(mpm->prho, 0, gv->NX, gv->NY, 7, 0, "prho");
    debug_check_matrix(mpm->pc11, 0, gv->NX, gv->NY, 7, 0, "pc11");
    debug_check_matrix(mpm->pc33, 0, gv->NX, gv->NY, 7, 0, "pc33");
    debug_check_matrix(mpm->pc13, 0, gv->NX, gv->NY, 7, 0, "pc13");
#endif

    return;
}
