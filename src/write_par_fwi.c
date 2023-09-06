
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

/*----------------------------------------------------------------------*
 * Write summary of parameters to stdout                                *
 *----------------------------------------------------------------------*/

#include "fd.h"
#include "macros.h"
#include "logging.h"
#include <stdbool.h>

void write_par_fwi(GlobVar *gv, GlobVarInv *vinv)
{
    int l;
    //char buffer[STRING_SIZE]="", buf[STRING_SIZE];

    log_info("----------------  FWI specific parameters  ------------------\n");

    log_info("Maximum number of iterations: %d\n", vinv->ITERMAX);
    log_info("Location of the measured seismograms: %s\n", vinv->DATA_DIR);
    switch (vinv->ADJOINT_TYPE) {
      case 1:
          log_info("Inversion of x and y component.\n");
          break;
      case 2:
          log_info("Inversion of y component only.\n");
          break;
      case 3:
          log_info("Inversion of x component only.\n");
          break;
      case 4:
          log_info("Inversion of pressure component.\n");
          break;
      default:
          log_fatal("Sorry, incorrect specification of parameter ADJOINT_TYPE!\n");
    }

    if (vinv->VELOCITY == 1) {
        log_info("Minimization of misfit in velocity seismograms.\n");
    } else {
        log_info("Minimization of misfit in displacement seismograms.\n");
    }

    log_info("Shots used for step length estimation (in total %d testshots are used):\n", vinv->NO_OF_TESTSHOTS);
    log_info("\tFirst testshot: %d \n", vinv->TESTSHOT_START);
    log_info("\tLast testshot: %d \n", vinv->TESTSHOT_END);
    log_info("\tTestshot increment: %d \n", vinv->TESTSHOT_INCR);

    log_info("Log file for misfit in each iteration step: %s\n", vinv->MISFIT_LOG_FILE);

    log_info("File name for output of inverted models: %s\n", vinv->INV_MODELFILE);
    log_info("Every %d iteration(s) a model will be output starting at iteration %d.\n", vinv->NF, vinv->NFSTART);
    log_info("File name for output of jacobian: %s\n", vinv->JACOBIAN);
    log_info("Every %d iteration(s) a jacobian matrix will be output starting at iteration %d.\n", vinv->NF_JAC,
             vinv->NFSTART_JAC);

    log_info("Vp is inverted from iteration step %d on.\n", vinv->INV_VP_ITER);
    if (gv->WEQ > 2 && gv->WEQ < 9) {
        log_info("Vs is inverted from iteration step %d on.\n", vinv->INV_VS_ITER);
    }
    log_info("Density is inverted from iteration step %d on.\n", vinv->INV_RHO_ITER);

    if (vinv->VP_VS_RATIO < 1) {
        log_info("Minimum Vp/Vs-ratio is <1 which means that it is disregarded.\n");
    } else {
        log_info("Minimum Vp/Vs-ratio is set to %4.2f.\n", vinv->VP_VS_RATIO);
    }

    if (vinv->S == 1) {
        log_info("Limited update of Vp in reference to the starting model is set to %4.2f %%.\n", vinv->S_VP);
        if (gv->WEQ > 2 && gv->WEQ < 9) {
            log_info("Limited update of Vs in reference to the starting model is set to %4.2f %%.\n", vinv->S_VS);
        }
        log_info("Limited update of Rho in reference to the starting model is set to %4.2f %%.\n", vinv->S_RHO);
    }

    log_info("------------------- Gradient tapering -----------------------\n");
    if (vinv->SWS_TAPER_GRAD_VERT == 1) {
        log_info("Vertical taper applied.\n");
    } else {
        log_info("No vertical taper applied.\n");
    }
    if (vinv->SWS_TAPER_GRAD_HOR == 1) {
        log_info("Horizontal taper applied.\n");
    } else {
        log_info("No horizontal taper applied.\n");
    }
    if (vinv->SWS_TAPER_GRAD_VERT == 1 || vinv->SWS_TAPER_GRAD_HOR == 1) {
        log_info
            ("Corner points defining slope of vertical and/or horizontal filter: GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d\n",
             vinv->GRADT1, vinv->GRADT2, vinv->GRADT3, vinv->GRADT4);
    }

    if (vinv->SWS_TAPER_GRAD_SOURCES == 1) {
        log_info("Taper around source applied.\n");
    }
    if (vinv->SWS_TAPER_CIRCULAR_PER_SHOT == 1) {
        log_info("Taper around source for each shot applied.\n");
    }
    if (vinv->SWS_TAPER_GRAD_SOURCES == 1 || vinv->SWS_TAPER_CIRCULAR_PER_SHOT == 1) {
        if (vinv->SRTSHAPE == 1) {
            log_info("Shape of taper is defined by an error function with a radius of %f.)\n", vinv->SRTRADIUS);
        } else if (vinv->SRTSHAPE == 2) {
            log_info("Shape of taper is defined by a log function with a radius of %f.)\n", vinv->SRTRADIUS);
        } else {
            log_fatal("Sorry, incorrect specification of parameter SRTSHAPE!\n");
        }
        if (vinv->FILTSIZE > 0) {
            log_info("An additional zero-taper of size 2*%d+1 x 2*%d+1 grid points is applied.\n", vinv->FILTSIZE,
                     vinv->FILTSIZE);
        }
    } else {
        log_info("No taper around the sources applied.\n");
    }

    if (vinv->SWS_TAPER_FILE == 1) {
        log_info("Taper files taper.bin, taper_u.bin taper_rho.bin are read in and applied to the gradients.\n");
    } else {
        log_info("No taper files are applied to the summed gradients.\n");
    }

    if (vinv->SWS_TAPER_FILE_PER_SHOT == 1) {
        log_info("Taper files for single shots are read in and applied to the gradients.\n");
        log_info("\tFile for vp gradients: %s.vp\n", vinv->TAPER_FILE_NAME);
        if (gv->WEQ > 2 && gv->WEQ < 9) {
            log_info("\tFile for vs gradients: %s.vs\n", vinv->TAPER_FILE_NAME);
        }
        log_info("\tFile for rho gradients: %s.rho\n", vinv->TAPER_FILE_NAME);
    } else {
        log_info("No taper files are applied to the gradients before summation.\n");
    }

    if (vinv->SPATFILTER == 1) {
        log_info("Smoothing (spatial filtering) of the gradients applied.\n");
        log_info("Size of smoothing filter (SPAT_FILT_SIZE): %d", vinv->SPAT_FILT_SIZE);
        log_info("Additional parameter of smoothing (SPAT_FILT_1): %d\n", vinv->SPAT_FILT_1);
        log_info("Additional parameter of smoothing (SPAT_FILT_ITER): %d)\n", vinv->SPAT_FILT_ITER);
    } else {
        log_info("No smoothing (spatial filtering) of the gradients applied.\n");
    }

    log_info("-------- Gradient smoothing with 2D-Gaussian filter ---------\n");
    if (vinv->GRAD_FILTER == 1) {
        if (vinv->GRAD_FILT_WAVELENGTH == 0)
            log_info("Gradients are filtered with a filter size of %d)\n", vinv->FILT_SIZE_GRAD);
        if (vinv->GRAD_FILT_WAVELENGTH == 1)
            log_info
                ("FILT_SIZE_GRAD is ignored. Gradients are filtered with a wavelength dependent filter size. Weighting factor A = %4.2f \n",
                 vinv->A);
    } else {
        log_info("Jacobians are not filtered.\n");
    }

    log_info("--------------- Limits of model parameters ------------------\n");
    log_info("VPLOWERLIM = %f \t VPUPPERLIM = %f \n", vinv->VPLOWERLIM, vinv->VPUPPERLIM);
    log_info("VSLOWERLIM = %f \t VSUPPERLIM = %f \n", vinv->VSLOWERLIM, vinv->VSUPPERLIM);
    log_info("RHOLOWERLIM = %f \t RHOUPPERLIM = %f \n", vinv->RHOLOWERLIM, vinv->RHOUPPERLIM);

    log_info("---- Calculation of diagonal elements of approx. Hessian ----\n");
    switch (vinv->GRAD_METHOD) {
      case 1:
          log_info("Applied Gradient method: PCG\n", vinv->GRAD_METHOD);
          break;
      case 2:
          log_info("Applied Gradient method: LBFGS\n", vinv->GRAD_METHOD);
          break;
      case 0:
          break;                /* only forward modeling is applied */
      default:
          log_fatal("Sorry, incorrect value for GRAD_METHOD!\n");
    }

    log_info("--------------------- Model smoothing -----------------------\n");
    if (vinv->MODEL_FILTER == 1) {
        log_info("Vp and Vs models are filtered after each iteration step (MODEL_FILTER=%d).\n", vinv->MODEL_FILTER);
        log_info("Used filter size: %d)\n", vinv->FILT_SIZE);
    } else {
        log_info("Vp and Vs models are not filtered after each iteration step (MODEL_FILTER=%d).\n",
                 vinv->MODEL_FILTER);
    }

    log_info("---------- Inversion of the source time function ------------\n");
    if (vinv->INV_STF == 1) {
        log_info("Source time function will be inverted (INV_STF=%d).\n", vinv->INV_STF);
        log_info("Parameter for STF: PARA=%s, N_STF=%d, N_STF_START=%d.\n", vinv->PARA, vinv->N_STF, vinv->N_STF_START);
    } else {
        log_info("No inversion of the source time function (INV_STF=%d).\n", vinv->INV_STF);
    }

    log_info("--------------------- Trace kill STF ------------------------\n");
    if (vinv->TRKILL_STF) {
        log_info("Trace kill for STF is applied (TRKILL_STF=%d). \n", vinv->TRKILL_STF);
        log_info("Reading trace kill STF matrix from file: %s \n\n", vinv->TRKILL_FILE_STF);
    } else {
        log_info("No trace kill STF is applied (TRKILL_STF=%d). \n", vinv->TRKILL_STF);
    }

    log_info("-------------------- Frequency filtering --------------------\n");
    if (vinv->TIME_FILT) {
        if (vinv->TIME_FILT == 1) {
            log_info("Time domain filtering is applied (TIME_FILT=%d).\n", vinv->TIME_FILT);
            log_info("Starting at frequencies of %.2f Hz.\n", vinv->F_LOW_PASS_START);
            log_info("Increasing the bandwidth up to %.2f Hz in steps of %.2f Hz.\n", vinv->F_LOW_PASS_END,
                     vinv->F_LOW_PASS_INCR);
        }
        if (vinv->TIME_FILT == 2) {
            log_info("Time domain filtering is applied (TIME_FILT=%d).\n", vinv->TIME_FILT);
            log_info("Frequencies will be read from file: %s\n", vinv->FREQ_FILE);
            log_info("Number of freqencies specified: %d \n", vinv->NFREQ);
            for (l = 1; l <= (vinv->NFREQ); l++) {
                log_info("Cut off frequency %d: %.2f Hz\n", l, vinv->F_LOW_PASS_EXT[l]);
            }
        }
        log_info("Order of lowpass filter is:\t%d\n", vinv->ORDER);
        if ((vinv->ORDER % 2) != 0) {
            log_fatal("Order of time domain filter must be an even number!\n");
        }
    } else {
        log_info("No time domain filtering is applied (TIME_FILT=%d).\n", vinv->TIME_FILT);
    }

    log_info("------------------------ Trace kill -------------------------\n");
    if (vinv->TRKILL) {
        log_info("Trace kill is applied (TRKILL=%d).\n", vinv->TRKILL);
        if (vinv->TRKILL_OFFSET) {
            log_info("Traces with offset between %f m and %f m are killed.\n", vinv->TRKILL_OFFSET_LOWER,
                     vinv->TRKILL_OFFSET_UPPER);
        } else {
            log_info("Reading trace kill matrix from file: %s \n", vinv->TRKILL_FILE);
        }
    } else {
        log_info("No trace kill is applied (TRKILL=%d).\n", vinv->TRKILL);
    }

    log_info("--------------- Time windowing and damping ------------------\n");
    if (vinv->TIMEWIN) {
        log_info("Time windowing and damping is applied (TIMEWIN=%d).\n", vinv->TIMEWIN);
        log_info("Reading picked times from files: %s\n", vinv->PICKS_FILE);
        log_info("Length of window after pick in s is: %f\n", vinv->TWLENGTH_PLUS);
        log_info("Length of window befor pick in s is: %f\n", vinv->TWLENGTH_MINUS);
        log_info("Gamma is: %f\n", vinv->GAMMA);
    } else {
        log_info("No time windowing and damping is applied (TIMEWIN=%d).\n", vinv->TIMEWIN);
    }

    log_info("------------------- Trace normalization ---------------------\n");
    if (vinv->NORMALIZE) {
        log_info("The measured and synthetic seismograms will be normalized\n", vinv->NORMALIZE);
        log_info(" before calculating the residuals (NORMALIZE=%d). \n");
    } else {
        log_info("No normalization of measured and synthetic seismograms (NORMALIZE=%d).\n", vinv->NORMALIZE);
    }

    log_info("------------------------- Workflow --------------------------\n");
    if (vinv->USE_WORKFLOW) {
        log_info("The workflow parameters are read from %s\n", vinv->FILE_WORKFLOW);
        log_info("Number of lines in WORKFLOW_FILE: %d\n", vinv->WORKFLOW_LINES);
        log_info("%s", vinv->WORKFLOW_HEADER);
        for (l = 1; l <= vinv->WORKFLOW_LINES; l++) {
            log_info("%3.0f %5.0f %5.0f %8.4f %5.0f %8.2f %8.2f %8.3f %8.3f %5.0f %9.4f %5.0f %10.2f\n",
                     (vinv->WORKFLOW)[l][1], (vinv->WORKFLOW)[l][2], (vinv->WORKFLOW)[l][3], (vinv->WORKFLOW)[l][4],
                     (vinv->WORKFLOW)[l][5], (vinv->WORKFLOW)[l][6], (vinv->WORKFLOW)[l][7], (vinv->WORKFLOW)[l][8],
                     (vinv->WORKFLOW)[l][9], (vinv->WORKFLOW)[l][10], (vinv->WORKFLOW)[l][11], (vinv->WORKFLOW)[l][12],
                     (vinv->WORKFLOW)[l][13], (vinv->WORKFLOW)[l][14]);
        }
    } else {
        log_info("No workflow provided.\n");
    }

    log_info("------------------- Gradient calculation --------------------\n");
    switch (vinv->LNORM) {
      case 1:
          log_info("L1 Norm is used.\n");
          break;
      case 2:
          log_info("L2 Norm is used.\n");
          break;
      case 3:
          log_info("Cauchy is used.\n");
          break;
      case 4:
          log_info("SECH is used.\n");
          break;
      case 5:
          log_info("Global correlation is used.\n");
          break;
      case 7:
          log_info("Normalized L2 Norm is used (each trace is normalized by its RMS value).\n");
          break;
      case 8:
          log_info("Enveloped-based Norm is used.\n");
          log_info("WATERLEVEL_LNORM8 is %e\n", vinv->WATERLEVEL_LNORM8);
          break;
      default:
          log_fatal("Sorry, incorrect value for LNORM!\n");
    }
    switch (vinv->DTINV) {
      case 0:
      case 1:
          log_info("Every time sample is used for the calculation of the gradients (DTINV).\n");
          break;
      case 2:
          log_info("Every %dnd time sample is used for the calculation of the gradients (DTINV).\n", vinv->DTINV);
          break;
      case 3:
          log_info("Every %drd time sample is used for the calculation of the gradients (DTINV).\n", vinv->DTINV);
          break;
      default:
          log_info("Every %dth time sample is used for the calculation of the gradients (DTINV).\n", vinv->DTINV);
    }

    log_info("------------------ Step length estimation -------------------\n");
    log_info("EPS_SCALE = %f\n", vinv->EPS_SCALE);
    log_info("STEPMAX = %d\n", vinv->STEPMAX);
    log_info("SCALEFAC = %f\n", vinv->SCALEFAC);

    log_info("-------------------- Termination criterium ------------------\n");
    log_info("Misfit change during the last two iterations is smaller than %f[%].\n", (vinv->PRO * 100.0));
    log_info("Minimum number of iteration per frequency (MIN_ITER) : %d \n", vinv->MIN_ITER);

    log_info("-------------------------------------------------------------\n");
    return;
}
