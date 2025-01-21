
/*------------------------------------------------------------------------
 * Copyright (C) 2016 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 3D-AWAIT. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 --------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * Calculation of L-BFGS update
 --------------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "macros.h"

//void lbfgs_reset(int iter, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS1, float ** y_LBFGS1, float * rho_LBFGS1);
//void lbfgs_reset(int iter, MemInv * minv, GlobVar *gv, GlobVarInv *vinv);
//void lbfgs_core(int iteration, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS, float ** y_LBFGS, float * rho_LBFGS,float *q_LBFGS,float *alpha_LBFGS,float *r_LBFGS);
//void lbfgs_core(int iter, float *q_LBFGS, float *alpha_LBFGS, float *r_LBFGS, MemInv * minv, GlobVar *gv, GlobVarInv *vinv);

void lbfgs(float **grad_vs, float **grad_rho, float **grad_vp, float Vs_avg, float rho_avg, float Vp_avg, int iter,
           MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{

    /* local */
    int w = 0;
    int i, j, l;
    float *q_LBFGS, *alpha_LBFGS, *r_LBFGS;
    //char grad[225];
    //FILE *FP_GRAD = NULL;

    /*---------------------*/
    /*      DEBUGGING      */

    /*---------------------*/

    /* if(!ACOUSTIC) {
     * sprintf(grad,"%s_grad1_vs_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_vs, grad);
     * }
     * 
     * sprintf(grad,"%s_grad1_rho_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_rho, grad);
     * 
     * if(WAVETYPE==1||WAVETYPE==3) {
     * sprintf(grad,"%s_grad1_vp_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_vp, grad);
     * } */

    /*---------------------*/
    /*      Experimental   */

    /*---------------------*/
    /* Scale the gradients with a constant factor */
    /* This is to avoid numerical instabilities if the gradients are to small (absolute value) */
    /* Do not use this feature, unless you have to */
    if (vinv->LBFGS_SCALE_GRADIENTS != 1) {
        if (gv->MPID == 0) {
            log_info("Scaling the gradients to ensure L-BFGS stability.\n");
            log_info("Scaling with factor %f\n", vinv->LBFGS_SCALE_GRADIENTS);
            log_info("This is an experimental feature.");
        }

        for (i = 1; i <= gv->NX; i++) {
            for (j = 1; j <= gv->NY; j++) {
                grad_vs[j][i] = grad_vs[j][i] * vinv->LBFGS_SCALE_GRADIENTS;
                if (LBFGS_NPAR > 1)
                    grad_rho[j][i] = grad_rho[j][i] * vinv->LBFGS_SCALE_GRADIENTS;
                if (LBFGS_NPAR > 2)
                    grad_vp[j][i] = grad_vp[j][i] * vinv->LBFGS_SCALE_GRADIENTS;
            }
        }
    }

    /*-------------------------------------------------*/
    /*      Init L-BFGS at iter==LBFGS_ITER_START      */

    /*-------------------------------------------------*/
    if (iter == vinv->LBFGS_ITER_START) {
        w = iter % vinv->N_LBFGS;
        if (w == 0)
            w = vinv->N_LBFGS;

        l = 0;
        for (i = 1; i <= gv->NX; i++) {
            for (j = 1; j <= gv->NY; j++) {
                l++;
                minv->y_LBFGS[w][l] = -grad_vs[j][i] * Vs_avg;  /* VS */
                if (LBFGS_NPAR > 1)
                    minv->y_LBFGS[w][l + gv->NX * gv->NY] = -grad_rho[j][i] * rho_avg;  /* RHO */
                if (LBFGS_NPAR > 2)
                    minv->y_LBFGS[w][l + 2 * gv->NX * gv->NY] = -grad_vp[j][i] * Vp_avg;    /* VP */
            }
        }
    }

    /*------------------------*/
    /*      Start L-BFGS      */

    /*------------------------*/
    if (iter > vinv->LBFGS_ITER_START) {

        alpha_LBFGS = vector(1, vinv->N_LBFGS);
        q_LBFGS = vector(1, LBFGS_NPAR * gv->NX * gv->NY);
        r_LBFGS = vector(1, LBFGS_NPAR * gv->NX * gv->NY);

        w = (iter - 1) % vinv->N_LBFGS;
        if (w == 0)
            w = vinv->N_LBFGS;

        /* Debugging */
        /*if(!ACOUSTIC) {
         * sprintf(grad,"%s_y_LBFGS_vs_it%d.bin.%i.%i",GRADIENT,iteration,POS[1],POS[2]);
         * FP_GRAD=fopen(grad,"wb");
         * } */

        if (gv->MPID == 0) {
            log_info("--------------- L-BFGS ---------------\n");
            log_info("Start calculation L-BFGS update.\n");
            log_info("At Iteration %i in L-BFGS vector %i\n", iter, w);
        }

        l = 0;
        for (i = 1; i <= gv->NX; i++) {
            for (j = 1; j <= gv->NY; j++) {
                l++;
                /* VS */
                //if(!ACOUSTIC){
                minv->y_LBFGS[w][l] += grad_vs[j][i] * Vs_avg;  /* add grad(i) to build grad(i)-grad(i-1) */
                q_LBFGS[l] = grad_vs[j][i] * Vs_avg;    /* Normalisation */
                //}
                /* RHO */
                if (LBFGS_NPAR > 1) {
                    minv->y_LBFGS[w][l + gv->NY * gv->NX] += grad_rho[j][i] * rho_avg;  /* add grad(i) to build grad(i)-grad(i-1) */
                    q_LBFGS[l + gv->NY * gv->NX] = grad_rho[j][i] * rho_avg;    /* Normalisation */
                }
                /* VP */
                if (LBFGS_NPAR > 2) {
                    minv->y_LBFGS[w][l + 2 * gv->NY * gv->NX] += grad_vp[j][i] * Vp_avg;    /* add grad(i) to build grad(i)-grad(i-1) */
                    q_LBFGS[l + 2 * gv->NY * gv->NX] = grad_vp[j][i] * Vp_avg;  /* Normalisation */
                }

                /* Debugging */
                /*if(!ACOUSTIC) {
                 * fwrite(&y_LBFGS[w][l],sizeof(float),1,FP_GRAD);
                 * } */
            }
        }

        /*---------------------*/
        /*      DEBUGGING      */

        /*---------------------*/
        /*if(!ACOUSTIC) {
         * fclose(FP_GRAD);
         * MPI_Barrier(MPI_COMM_WORLD);
         * sprintf(grad,"%s_y_LBFGS_vs_it%d.bin",GRADIENT,iteration);
         * if (gv->MPID==0) mergemod(grad,3);
         * MPI_Barrier(MPI_COMM_WORLD);
         * sprintf(grad,"%s_y_LBFGS_vs_it%d.bin.%i.%i",GRADIENT,iteration,POS[1],POS[2]);
         * remove(grad);
         * } */

        /*----------------------------------*/
        /*      call L-BFGS Algorithm       */

        /*----------------------------------*/
        lbfgs_core(iter, q_LBFGS, alpha_LBFGS, r_LBFGS, minv, gv, vinv);

        /*-------------------------------------------------------------*/
        /* Save model pertubation and save gradient for next iteration */

        /*-------------------------------------------------------------*/
        w = iter % vinv->N_LBFGS;
        if (w == 0)
            w = vinv->N_LBFGS;

        l = 0;
        for (i = 1; i <= gv->NX; i++) {
            for (j = 1; j <= gv->NY; j++) {
                l++;
                /* VS */
                //if(!ACOUSTIC) {
                minv->y_LBFGS[w][l] = -grad_vs[j][i] * Vs_avg;  /* add -grad(i-1) to build grad(i)-grad(i-1) */
                grad_vs[j][i] = r_LBFGS[l] * Vs_avg;    /* Denormalization */
                //}
                /* RHO */
                //if(LBFGS_NPAR>1) {
                minv->y_LBFGS[w][l + gv->NY * gv->NX] = -grad_rho[j][i] * rho_avg;  /* add -grad(i-1) to build grad(i)-grad(i-1) */
                grad_rho[j][i] = r_LBFGS[l + gv->NY * gv->NX] * rho_avg;    /* Denormalization */
                //}
                /* VP */
                //if(LBFGS_NPAR>2) {
                minv->y_LBFGS[w][l + 2 * gv->NY * gv->NX] = -grad_vp[j][i] * Vp_avg;    /* add -grad(i-1) to build grad(i)-grad(i-1) */
                grad_vp[j][i] = r_LBFGS[l + 2 * gv->NY * gv->NX] * Vp_avg;  /* Denormalization */
                //}
            }
        }

        free_vector(r_LBFGS, 1, LBFGS_NPAR * gv->NX * gv->NY);
        free_vector(q_LBFGS, 1, LBFGS_NPAR * gv->NX * gv->NY);
        free_vector(alpha_LBFGS, 1, vinv->N_LBFGS);

    }

    /*---------------------*/
    /*      DEBUGGING      */

    /*---------------------*/
    /*if(!ACOUSTIC){
     * sprintf(grad,"%s_grad2_vs_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_vs, grad);
     * }
     * 
     * sprintf(grad,"%s_grad2_rho_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_rho, grad);
     * 
     * if(WAVETYPE==1||WAVETYPE==3) {
     * sprintf(grad,"%s_grad2_vp_it%d",GRADIENT,iteration);
     * write_matrix_disk(grad_vp, grad);
     * } */
}

void lbfgs_core(int iter, float *q_LBFGS, float *alpha_LBFGS, float *r_LBFGS, MemInv * minv, GlobVar *gv,
                GlobVarInv *vinv)
{

    float beta_LBFGS = 0.0;
    float dum1 = 0.0, dum2 = 0.0, buf1 = 0.0, buf2 = 0.0;
    float h0;
    int v = 0, w = 0, l = 0;
    int lmax = gv->NX * gv->NY * LBFGS_NPAR;
    int m = iter - vinv->N_LBFGS;
    if (m < 1)
        m = 1;

    /*----------------------------------*/
    /* calculate H0 and rho_LBFGS       */

    /*----------------------------------*/
    w = (iter - 1) % vinv->N_LBFGS;
    if (w == 0)
        w = vinv->N_LBFGS;
    dum1 = 0.0;
    dum2 = 0.0;

    for (l = 1; l <= lmax; l++) {
        dum1 += minv->y_LBFGS[w][l] * minv->s_LBFGS[w][l];
        dum2 += minv->y_LBFGS[w][l] * minv->y_LBFGS[w][l];
    }

    buf1 = 0.0;
    buf2 = 0.0;
    MPI_Allreduce(&dum1, &buf1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&dum2, &buf2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    minv->rho_LBFGS[w] = 1 / buf1;

    h0 = buf1 / buf2;

    /* give output so stdout */
    /*if(VERBOSE || 1) {
     * for(w=1;w<=vinv->N_LBFGS;w++) {
     * fprintf(FP,"\n rho_LBFGS(%2d)=%e",w,minv->rho_LBFGS[w]);
     * }
     * fprintf(FP,"\n h0=%e\n",h0);
     * } */

    /*----------------------------------*/
    /*       L-BFGS loop 1              */

    /*----------------------------------*/

    for (v = iter - 1; v >= m; v--) {

        w = v % vinv->N_LBFGS;
        if (w == 0)
            w = vinv->N_LBFGS;

        alpha_LBFGS[w] = 0.0;
        for (l = 1; l <= LBFGS_NPAR; l++) {
            alpha_LBFGS[w] += minv->rho_LBFGS[w] * minv->s_LBFGS[w][l] * q_LBFGS[l];
        }

        buf1 = 0.0;
        dum2 = alpha_LBFGS[w];
        MPI_Allreduce(&dum2, &buf1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        alpha_LBFGS[w] = buf1;

        for (l = 1; l <= LBFGS_NPAR; l++) {
            q_LBFGS[l] = q_LBFGS[l] - alpha_LBFGS[w] * minv->y_LBFGS[w][l];
        }
    }

    /*----------------------------------*/
    /*       Apply H0^-1                */

    /*----------------------------------*/
    for (l = 1; l <= LBFGS_NPAR; l++) {
        r_LBFGS[l] = h0 * q_LBFGS[l];
    }

    /*----------------------------------*/
    /*       L-BFGS loop 2              */

    /*----------------------------------*/
    for (v = m; v <= iter - 1; v++) {

        w = v % vinv->N_LBFGS;
        if (w == 0)
            w = vinv->N_LBFGS;

        beta_LBFGS = 0.0;

        for (l = 1; l <= LBFGS_NPAR; l++) {
            beta_LBFGS += minv->rho_LBFGS[w] * minv->y_LBFGS[w][l] * r_LBFGS[l];
        }

        buf1 = 0.0;
        buf2 = beta_LBFGS;
        MPI_Allreduce(&buf2, &buf1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        beta_LBFGS = buf1;

        for (l = 1; l <= LBFGS_NPAR; l++) {
            r_LBFGS[l] = r_LBFGS[l] + minv->s_LBFGS[w][l] * (alpha_LBFGS[w] - beta_LBFGS);
        }
    }

}

//void lbfgs_reset(int iter, int N_LBFGS, int LBFGS_NPAR,float ** s_LBFGS1, float ** y_LBFGS1, float * rho_LBFGS1);
void lbfgs_reset(int iter, MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{

    /* local variables */
    int l, m;

    if (gv->MPID == 0) {
        log_info("--------------- L-BFGS ---------------\n");
        if (iter > 1) {
            log_info("Reset L-BFGS at iteration %d.\n", iter);
        } else if (iter == 1) {
            log_info("L-BFGS will be used from iteration %d on.\n", iter + 1);
        }
    }

    for (l = 1; l <= vinv->N_LBFGS; l++) {
        for (m = 1; m <= (LBFGS_NPAR * gv->NX * gv->NY); m++) {
            minv->s_LBFGS[l][m] = 0.0;
            minv->y_LBFGS[l][m] = 0.0;
        }
        minv->rho_LBFGS[l] = 0.0;
    }
}
