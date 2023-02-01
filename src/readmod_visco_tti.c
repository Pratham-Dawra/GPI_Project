
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

/* --------------------------------------------------------------------------------------
 * Read viscoelastic TTI model properties (vp,vs,rho,eps,delta,theta,qp,qs) from files
 * ------------------------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "read_su.h"
#include <unistd.h>
#include <stdbool.h>

#ifdef EBUG
#include "su_struct.h"
#include "write_su.h"
#define NOUT 13
#endif

void readmod_visco_tti(MemModel * mpm, GlobVar * gv)
{
    float c11, c33, c13, c55, tau11, tau33, tau13, tau55, tau15, tau35;
    float *pts, sumc11, sumc13, sumc33, sumc55;
    float l1, l2, l12, l22, l14, l24, l13, l23;
    float a1, a3, a4, a5, a6, t;
    float c11t, c33t, c55t, c13t, c15t, c35t;

    int ii, jj;
    size_t ny;
    char filename[STRING_SIZE + 16];
    bool b_issu = false;

    const char *model[] = { "P-wave velocity model", "S-wave velocity model", "density model",
        "Epsilon model", "Delta model", "Theta model", "Quality factor Qp model", "Quality factor Qs model"
    };
    const char *suffix[] = { "vp", "vs", "rho", "epsilon", "delta", "theta", "qp", "qs" };
    const char *suffix_su[] = { "vp.su", "vs.su", "rho.su", "epsilon.su", "delta.su", "theta.su", "qp.su", "qs.su" };

    enum PARAMOD {
        P_VP = 0,
        P_VS,
        P_RHO,
        P_EPS,
        P_DEL,
        P_TET,
        P_QP,
        P_QS,
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

    /* vector for maxwellbodies */
    pts = vector(1, gv->L);
    for (int l = 1; l <= gv->L; l++) {
        pts[l] = 1.0 / (2.0 * PI * gv->FL[l]);
        mpm->peta[l] = gv->DT / pts[l];
    }

    float fc = 1.0 / gv->TS;
    log_infoc(0, "TTI: center source frequency of %5.2fHz applied for calculation of relaxed moduli.\n", fc);

    float ws = 2.0 * PI * fc;

#ifdef EBUG
    // set up output of all initial matrices for debug purpose
    float **debug_buffer = NULL;
    FILE *fp_debug[NOUT];
    char modfile[STRING_SIZE + 16];
    const char *debug_name[] = { "pc11", "rho", "pc33", "pc13", "pc55", "pc15", "pc35", 
				 "ptau11", "ptau33", "ptau13", "ptau55", "ptau15", "ptau35" };

    if (0 == gv->MPID) {
	debug_buffer = (float **)malloc2d(NOUT, ny, sizeof(float));
	for (int i = 0; i < NOUT; ++i) {
	    sprintf(modfile, "%s.debug.%s.su", gv->MFILE, debug_name[i]);
	    log_debug("Opening debug file %s.\n", modfile);
	    fp_debug[i] = fopen(modfile, "wb");
	    if (!fp_debug[i])
		log_fatal("Could not open debug file %s for writing.\n", modfile);
	}
    }
#endif

    /* loop over global grid */
    for (int i = 1; i <= gv->NXG; i++) {
        if (b_issu) {
            su_read_trace(fp[P_VP], 0, (unsigned short)ny, false, NULL, &(para[P_VP][0]));
            su_read_trace(fp[P_VS], 0, (unsigned short)ny, false, NULL, &(para[P_VS][0]));
            su_read_trace(fp[P_RHO], 0, (unsigned short)ny, false, NULL, &(para[P_RHO][0]));
            su_read_trace(fp[P_EPS], 0, (unsigned short)ny, false, NULL, &(para[P_EPS][0]));
            su_read_trace(fp[P_DEL], 0, (unsigned short)ny, false, NULL, &(para[P_DEL][0]));
            su_read_trace(fp[P_TET], 0, (unsigned short)ny, false, NULL, &(para[P_TET][0]));
            su_read_trace(fp[P_QP], 0, (unsigned short)ny, false, NULL, &(para[P_QP][0]));
            su_read_trace(fp[P_QS], 0, (unsigned short)ny, false, NULL, &(para[P_QS][0]));
        } else {
            fread(&(para[P_VP][0]), sizeof(float), ny, fp[P_VP]);
            fread(&(para[P_VS][0]), sizeof(float), ny, fp[P_VS]);
            fread(&(para[P_RHO][0]), sizeof(float), ny, fp[P_RHO]);
            fread(&(para[P_EPS][0]), sizeof(float), ny, fp[P_EPS]);
            fread(&(para[P_DEL][0]), sizeof(float), ny, fp[P_DEL]);
            fread(&(para[P_TET][0]), sizeof(float), ny, fp[P_TET]);
            fread(&(para[P_QP][0]), sizeof(float), ny, fp[P_QP]);
            fread(&(para[P_QS][0]), sizeof(float), ny, fp[P_QS]);
        }
        for (int j = 1; j <= gv->NYG; j++) {
            /* calculation of required elastic constants:
             * c11, c33, c13, c55, tau11, tau33, tau55 */
            c33 = para[P_RHO][j - 1] * para[P_VP][j - 1] * para[P_VP][j - 1];
            c55 = para[P_RHO][j - 1] * para[P_VS][j - 1] * para[P_VS][j - 1];
            c11 = c33 * (2.0 * para[P_EPS][j - 1] + 1.0);
            c13 = sqrt((2.0 * para[P_DEL][j - 1] * c33 * (c33 - c55)) + ((c33 - c55) * (c33 - c55))) - c55;
            tau11 = tau33 = 2.0 / (para[P_QP][j - 1] * gv->L);
            tau55 = 2.0 / (para[P_QS][j - 1] * gv->L);
            sumc11 = 0.0f;
            sumc33 = 0.0f;
            sumc55 = 0.0f;
            for (int l = 1; l <= gv->L; l++) {
                sumc11 += ((ws * ws * pts[l] * pts[l] * tau11) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumc33 += ((ws * ws * pts[l] * pts[l] * tau33) / (1.0 + ws * ws * pts[l] * pts[l]));
                sumc55 += ((ws * ws * pts[l] * pts[l] * tau55) / (1.0 + ws * ws * pts[l] * pts[l]));
            }

            /* relaxed moduli */
            c11 = c11 / (1.0 + sumc11);
            c33 = c33 / (1.0 + sumc33);
            c55 = c55 / (1.0 + sumc55);

            /* isotropic attenuation */
            tau13 = tau15 = tau35 = (c11 * gv->L * tau11 - 2.0 * c55 * gv->L * tau55) / (c11 - 2.0 * c55);

            sumc13 = 0.0f;
            for (int l = 1; l <= gv->L; l++) {
                sumc13 = sumc13 + ((ws * ws * pts[l] * pts[l] * tau13) / (1.0 + ws * ws * pts[l] * pts[l]));
            }
            c13 = c13 / (1.0 + sumc13);

            /* Bond transformation (Oh et al, 2020, GJI, doi: 10.1093/gji/ggaa295 */
            t = para[P_TET][j - 1] * PI / 180.0;
            l1 = cos(t);
            l2 = sin(t);
            l12 = l1 * l1;
            l22 = l2 * l2;
            l14 = l12 * l12;
            l24 = l22 * l22;
            l13 = l1 * l12;
            l23 = l2 * l22;
            a1 = 2.0 * c13 + 4.0 * c55;
            a3 = c11 + c33 - 4.0 * c55;
            a4 = c11 + c33 - 2.0 * c13;
            a5 = c13 - c11 + 2.0 * c55;
            a6 = c13 - c33 + 2.0 * c55;
            c11t = c11 * l14 + c33 * l24 + a1 * l12 * l22;
            c33t = c11 * l24 + c33 * l14 + a1 * l12 * l22;
            c13t = a3 * l12 * l22 + c13 * (l14 + l24);
            c55t = a4 * l12 * l22 + c55 * (l12 - l22) * (l12 - l22);
            c15t = a5 * l13 * l2 - a6 * l1 * l23;
            c35t = a5 * l23 * l1 - a6 * l2 * l13;

#ifdef EBUG
	    if (0 == gv->MPID) {
		debug_buffer[0][j-1] = c11t;
		debug_buffer[1][j-1] = para[P_RHO][j - 1];
		debug_buffer[2][j-1] = c33t;
		debug_buffer[3][j-1] = c13t;
		debug_buffer[4][j-1] = c55t;
		debug_buffer[5][j-1] = c15t;
		debug_buffer[6][j-1] = c35t;
		debug_buffer[7][j-1] = tau11;
		debug_buffer[8][j-1] = tau33;
		debug_buffer[9][j-1] = tau13;
		debug_buffer[10][j-1] = tau55;
		debug_buffer[11][j-1] = tau15;
		debug_buffer[12][j-1] = tau35;
	    }
#endif

            /* only the PE which belongs to the current global gridpoint 
             * is saving model parameters in his local arrays */
            if ((gv->POS[1] == ((i - 1) / gv->NX)) && (gv->POS[2] == ((j - 1) / gv->NY))) {
                ii = i - gv->POS[1] * gv->NX;
                jj = j - gv->POS[2] * gv->NY;
                mpm->pc11[jj][ii] = c11t;
                mpm->prho[jj][ii] = para[P_RHO][j - 1];
                mpm->pc33[jj][ii] = c33t;
                mpm->pc13[jj][ii] = c13t;
                mpm->pc55[jj][ii] = c55t;
                mpm->pc15[jj][ii] = c15t;
                mpm->pc35[jj][ii] = c35t;
                mpm->ptau11[jj][ii] = tau11;
                mpm->ptau33[jj][ii] = tau33;
                mpm->ptau13[jj][ii] = tau13;
                mpm->ptau55[jj][ii] = tau55;
                mpm->ptau15[jj][ii] = tau15;
                mpm->ptau35[jj][ii] = tau35;
            }
        }

#ifdef EBUG
	// output matrix column (=trace)
	if (0 == gv->MPID) {
	    int ierr = 0;
	    SUhead head;
	    init_SUhead(&head);
	    head.ns =  (unsigned short)gv->NYG;
	    head.dt =  (int)gv->DH;
	    head.tracl = head.tracr = i;
	    head.ntr = gv->NXG;
	    head.d1 = head.d2 = gv->DH;	    
	    for (int k = 0; k < NOUT; ++k) {
		ierr = su_write_trace(fp_debug[k], &head, &(debug_buffer[k][0]));
		if (0 != ierr) {
		    sprintf(modfile, "%s.debug.%s.su", gv->MFILE, debug_name[k]);   
		    log_warn("Error while writing SU debug file %s.\n", modfile);
		}
	    }
	}
#endif

    }

    free_vector(pts, 1, gv->L);
    if (para)
        free(para);

#ifdef EBUG
    if (0 == gv->MPID) {
	if (debug_buffer)
	    free(debug_buffer);
	for (int i = 0; i < NOUT; ++i) {
	    fclose(fp_debug[i]);
	}
    }
#endif

    for (int i = 0; i < NPARA; ++i) {
        fclose(fp[i]);
    }

    return;
}
