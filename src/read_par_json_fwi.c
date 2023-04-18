
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

/*------------------------------------------------------------------------
 *   program SOFI2D, reading FWI input-parameters from input-file or stdin
 *   that are formatted according to the json standard
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"
#include "logging.h"

//char ** varname_list,** value_list;

int read_par_json_fwi(int number_readobjects, char **varname_list, char **value_list, int *used_list, int fserr,
                      GlobVar *gv, GlobVarInv *vinv)
{

    /* extract variables form object list */

    /*=================================
    section source time function inversion (STF)
    =================================*/
    log_warn("The following default values are set:\n");
    log_warn("=====================================\n");

    /* Definition of inversion for source time function */
    if (get_int_from_objectlist("INV_STF", number_readobjects, &(vinv->INV_STF), varname_list, value_list, used_list)) {
        vinv->INV_STF = 0;
        log_warn("Variable INV_STF is set to default value %d.\n", vinv->INV_STF);

    } else {
        if (vinv->INV_STF == 1) {
            if (get_string_from_objectlist
                ("PARA", number_readobjects, (vinv->PARA), varname_list, value_list, used_list))
                log_fatal("Variable PARA could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("N_STF", number_readobjects, &(vinv->N_STF), varname_list, value_list, used_list))
                log_fatal("Variable N_STF could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("N_STF_START", number_readobjects, &(vinv->N_STF_START), varname_list, value_list, used_list))
                log_fatal("Variable N_STF_START could not be retrieved from the json input file!");
        }
    }

    /*=================================
    section inversion parameters
    =================================*/

    /* General inversion parameters */
    if (get_int_from_objectlist("ITERMAX", number_readobjects, &(vinv->ITERMAX), varname_list, value_list, used_list))
        log_fatal("Variable ITERMAX could not be retrieved from the json input file!");
    if (get_string_from_objectlist
        ("DATA_DIR", number_readobjects, (vinv->DATA_DIR), varname_list, value_list, used_list))
        log_fatal("Variable DATA_DIR could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("ADJOINT_TYPE", number_readobjects, &(vinv->ADJOINT_TYPE), varname_list, value_list, used_list))
        log_fatal("Variable ADJOINT_TYPE could not be retrieved from the json input file!");
    if (get_string_from_objectlist
        ("MISFIT_LOG_FILE", number_readobjects, (vinv->MISFIT_LOG_FILE), varname_list, value_list, used_list)) {
        strcpy(vinv->MISFIT_LOG_FILE, "L2_LOG.dat");
        log_warn("Variable MISFIT_LOG_FILE is set to default value %s.\n", vinv->MISFIT_LOG_FILE);
    }
    if (get_int_from_objectlist("VELOCITY", number_readobjects, &(vinv->VELOCITY), varname_list, value_list, used_list)) {
        vinv->VELOCITY = 0;
        log_warn("Variable VELOCITY is set to default value %d.\n", vinv->VELOCITY);
    }
    if (get_int_from_objectlist
        ("USE_WORKFLOW", number_readobjects, &(vinv->USE_WORKFLOW), varname_list, value_list, used_list)) {
        vinv->USE_WORKFLOW = 0;
    }
    if (vinv->USE_WORKFLOW != 0) {
        if (get_string_from_objectlist
            ("FILE_WORKFLOW", number_readobjects, (vinv->FILE_WORKFLOW), varname_list, value_list, used_list))
            log_fatal("Variable FILE_WORKFLOW could not be retrieved from the json input file!");
    }

    if (get_int_from_objectlist("EPRECOND", number_readobjects, &(vinv->EPRECOND), varname_list, value_list, used_list)) {
        vinv->EPRECOND = 0;
        log_warn("Variable EPRECOND is set to default value %d.\n", vinv->EPRECOND);
    } else {
        if (get_int_from_objectlist
            ("EPRECOND_ITER", number_readobjects, &(vinv->EPRECOND_ITER), varname_list, value_list, used_list)) {
            vinv->EPRECOND_ITER = 0;
            log_warn("Variable EPRECOND_ITER is set to default value %d.\n", vinv->EPRECOND_ITER);
        }
        if (get_int_from_objectlist
            ("EPRECOND_PER_SHOT", number_readobjects, &(vinv->EPRECOND_PER_SHOT), varname_list, value_list,
             used_list)) {
            vinv->EPRECOND_PER_SHOT = 0;
            log_warn("Variable EPRECOND_PER_SHOT is set to default value %d.\n", vinv->EPRECOND_PER_SHOT);
            if (vinv->EPRECOND_ITER != 0) {
                vinv->EPRECOND_ITER = 0;
                log_warn(" EPRECOND_PER_SHOT and EPRECOND_ITER>0 not supported.\n");
                log_warn(" EPRECOND_ITER is set to EPRECOND_ITER=%d.\n", vinv->EPRECOND_ITER);
            }
        }
        if (get_float_from_objectlist
            ("EPSILON_WE", number_readobjects, &(vinv->EPSILON_WE), varname_list, value_list, used_list))
            log_fatal("Variable EPSILON_WE could not be retrieved from the json input file!");
    }
    vinv->EPRECOND_MAX = vinv->EPRECOND;

    if (get_int_from_objectlist
        ("TESTSHOT_START", number_readobjects, &(vinv->TESTSHOT_START), varname_list, value_list, used_list))
        log_fatal("Variable TESTSHOT_START could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("TESTSHOT_END", number_readobjects, &(vinv->TESTSHOT_END), varname_list, value_list, used_list))
        log_fatal("Variable TESTSHOT_START could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("TESTSHOT_INCR", number_readobjects, &(vinv->TESTSHOT_INCR), varname_list, value_list, used_list))
        log_fatal("Variable TESTSHOT_INCR could not be retrieved from the json input file!");
    /* calculation of number of testsshots */
    vinv->NO_OF_TESTSHOTS = (vinv->TESTSHOT_END - vinv->TESTSHOT_START) / vinv->TESTSHOT_INCR + 1;

    /* Definition of gradient taper geometry */
    if (get_int_from_objectlist
        ("SWS_TAPER_GRAD_VERT", number_readobjects, &(vinv->SWS_TAPER_GRAD_VERT), varname_list, value_list,
         used_list)) {
        vinv->SWS_TAPER_GRAD_VERT = 0;
        log_warn("Variable SWS_TAPER_GRAD_VERT is set to default value %d.\n", vinv->SWS_TAPER_GRAD_VERT);
    }
    if (get_int_from_objectlist
        ("SWS_TAPER_GRAD_HOR", number_readobjects, &(vinv->SWS_TAPER_GRAD_HOR), varname_list, value_list, used_list)) {
        vinv->SWS_TAPER_GRAD_HOR = 0;
        log_warn("Variable SWS_TAPER_GRAD_HOR is set to default value %d.\n", vinv->SWS_TAPER_GRAD_HOR);
    }
    if ((vinv->SWS_TAPER_GRAD_VERT == 1) || (vinv->SWS_TAPER_GRAD_HOR == 1)) {
        if (get_int_from_objectlist("GRADT1", number_readobjects, &(vinv->GRADT1), varname_list, value_list, used_list))
            log_fatal("Variable GRADT1 could not be retrieved from the json input file!");
        if (get_int_from_objectlist("GRADT2", number_readobjects, &(vinv->GRADT2), varname_list, value_list, used_list))
            log_fatal("Variable GRADT2 could not be retrieved from the json input file!");
        if (get_int_from_objectlist("GRADT3", number_readobjects, &(vinv->GRADT3), varname_list, value_list, used_list))
            log_fatal("Variable GRADT3 could not be retrieved from the json input file!");
        if (get_int_from_objectlist("GRADT4", number_readobjects, &(vinv->GRADT4), varname_list, value_list, used_list))
            log_fatal("Variable GRADT4 could not be retrieved from the json input file!");
    }
    if (get_int_from_objectlist
        ("SWS_TAPER_GRAD_SOURCES", number_readobjects, &(vinv->SWS_TAPER_GRAD_SOURCES), varname_list, value_list,
         used_list)) {
        vinv->SWS_TAPER_GRAD_SOURCES = 0;
        log_warn("Variable SWS_TAPER_GRAD_SOURCES is set to default value %d.\n", vinv->SWS_TAPER_GRAD_SOURCES);
    }
    if (get_int_from_objectlist
        ("SWS_TAPER_CIRCULAR_PER_SHOT", number_readobjects, &(vinv->SWS_TAPER_CIRCULAR_PER_SHOT), varname_list,
         value_list, used_list)) {
        vinv->SWS_TAPER_CIRCULAR_PER_SHOT = 0;
        log_warn("Variable SWS_TAPER_CIRCULAR_PER_SHOT is set to default value %d.\n",
                 vinv->SWS_TAPER_CIRCULAR_PER_SHOT);
    }
    if ((vinv->SWS_TAPER_GRAD_SOURCES == 1) || (vinv->SWS_TAPER_CIRCULAR_PER_SHOT == 1)) {
        if (get_int_from_objectlist
            ("SRTSHAPE", number_readobjects, &(vinv->SRTSHAPE), varname_list, value_list, used_list))
            log_fatal("Variable SRTSHAPE could not be retrieved from the json input file!");
        if (get_float_from_objectlist
            ("SRTRADIUS", number_readobjects, &(vinv->SRTRADIUS), varname_list, value_list, used_list))
            log_fatal("Variable SRTRADIUS could not be retrieved from the json input file!");
        if (get_int_from_objectlist
            ("FILTSIZE", number_readobjects, &(vinv->FILTSIZE), varname_list, value_list, used_list))
            log_fatal("Variable FILTSIZE could not be retrieved from the json input file!");
    }
    if (get_int_from_objectlist
        ("SWS_TAPER_FILE", number_readobjects, &(vinv->SWS_TAPER_FILE), varname_list, value_list, used_list)) {
        vinv->SWS_TAPER_FILE = 0;
        log_warn("Variable SWS_TAPER_FILE is set to default value %d.\n", vinv->SWS_TAPER_FILE);
    }
    if (vinv->SWS_TAPER_FILE == 1) {
        if (get_string_from_objectlist
            ("TAPER_FILE_NAME", number_readobjects, (vinv->TAPER_FILE_NAME), varname_list, value_list, used_list))
            log_fatal("Variable TAPER_FILE_NAME could not be retrieved from the json input file!");
    }
    if (get_int_from_objectlist
        ("SWS_TAPER_FILE_PER_SHOT", number_readobjects, &(vinv->SWS_TAPER_FILE_PER_SHOT), varname_list, value_list,
         used_list)) {
        vinv->SWS_TAPER_FILE_PER_SHOT = 0;
        log_warn("Variable SWS_TAPER_FILE_PER_SHOT is set to default value %d.\n", vinv->SWS_TAPER_FILE_PER_SHOT);
    }
    if (vinv->SWS_TAPER_FILE_PER_SHOT == 1) {
        if (get_string_from_objectlist
            ("TAPER_FILE_NAME", number_readobjects, (vinv->TAPER_FILE_NAME), varname_list, value_list, used_list))
            log_fatal("Variable TAPER_FILE_NAME could not be retrieved from the json input file!");
    }

    /* Definition of smoothing (spatial filtering) of the gradients */
    if (get_int_from_objectlist
        ("SPATFILTER", number_readobjects, &(vinv->SPATFILTER), varname_list, value_list, used_list)) {
        vinv->SPATFILTER = 0;
        log_warn("Variable SPATFILTER is set to default value %d.\n", vinv->SPATFILTER);
    } else {
        if (vinv->SPATFILTER == 1) {
            if (get_int_from_objectlist
                ("SPAT_FILT_SIZE", number_readobjects, &(vinv->SPAT_FILT_SIZE), varname_list, value_list, used_list))
                log_fatal("Variable SPAT_FILT_SIZE could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("SPAT_FILT_1", number_readobjects, &(vinv->SPAT_FILT_1), varname_list, value_list, used_list))
                log_fatal("Variable SPAT_FILT_1 could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("SPAT_FILT_ITER", number_readobjects, &(vinv->SPAT_FILT_ITER), varname_list, value_list, used_list))
                log_fatal("Variable SPAT_FILT_ITER could not be retrieved from the json input file!");
        }
    }

    /* Definition of 2D-Gaussian filter of the gradients */
    if (get_int_from_objectlist
        ("GRAD_FILTER", number_readobjects, &(vinv->GRAD_FILTER), varname_list, value_list, used_list)) {
        vinv->GRAD_FILTER = 0;
        log_warn("Variable GRAD_FILTER is set to default value %d.\n", vinv->GRAD_FILTER);
    }
    if (get_int_from_objectlist
        ("FILT_SIZE_GRAD", number_readobjects, &(vinv->FILT_SIZE_GRAD), varname_list, value_list, used_list)) {
        vinv->FILT_SIZE_GRAD = 0;
        log_warn("Variable FILT_SIZE_GRAD is set to default value %d.\n", vinv->FILT_SIZE_GRAD);
    }
    if (get_int_from_objectlist
        ("GRAD_FILT_WAVELENGTH", number_readobjects, &(vinv->GRAD_FILT_WAVELENGTH), varname_list, value_list,
         used_list)) {
        vinv->GRAD_FILT_WAVELENGTH = 0;
        log_warn("Variable GRAD_FILT_WAVELENGTH is set to default value %d.\n", vinv->GRAD_FILT_WAVELENGTH);
    } else {
        if (vinv->GRAD_FILT_WAVELENGTH == 1) {
            if (get_float_from_objectlist("A", number_readobjects, &(vinv->A), varname_list, value_list, used_list))
                log_fatal("Variable A could not be retrieved from the json input file!");
        }
    }

    /* Output of inverted models */
    if (get_string_from_objectlist
        ("INV_MODELFILE", number_readobjects, (vinv->INV_MODELFILE), varname_list, value_list, used_list))
        log_fatal("Variable INV_MODELFILE could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NFSTART", number_readobjects, &(vinv->NFSTART), varname_list, value_list, used_list))
        log_fatal("Variable NFSTART could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NF", number_readobjects, &(vinv->NF), varname_list, value_list, used_list))
        log_fatal("Variable NF could not be retrieved from the json input file!");

    /* Output of gradients */
    if (get_string_from_objectlist
        ("JACOBIAN", number_readobjects, (vinv->JACOBIAN), varname_list, value_list, used_list))
        log_fatal("Variable JACOBIAN could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("NFSTART_JAC", number_readobjects, &(vinv->NFSTART_JAC), varname_list, value_list, used_list))
        log_fatal("Variable NFSTART_JAC could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NF_JAC", number_readobjects, &(vinv->NF_JAC), varname_list, value_list, used_list))
        log_fatal("Variable NF_JAC could not be retrieved from the json input file!");

    /* Inversion for density */
    if (get_int_from_objectlist
        ("INV_RHO_ITER", number_readobjects, &(vinv->INV_RHO_ITER), varname_list, value_list, used_list)) {
        vinv->INV_RHO_ITER = 0;
        log_warn("Variable INV_RHO_ITER is set to default value %d.\n", vinv->INV_RHO_ITER);
    }
    /* Inversion for Vp */
    if (get_int_from_objectlist
        ("INV_VP_ITER", number_readobjects, &(vinv->INV_VP_ITER), varname_list, value_list, used_list)) {
        vinv->INV_VP_ITER = 0;
        log_warn("Variable INV_VP_ITER is set to default value %d.\n", vinv->INV_VP_ITER);
    }
    /* Inversion for Vs */
    if (get_int_from_objectlist
        ("INV_VS_ITER", number_readobjects, &(vinv->INV_VS_ITER), varname_list, value_list, used_list)) {
        vinv->INV_VS_ITER = 0;
        log_warn("Variable INV_VS_ITER is set to default value %d.\n", vinv->INV_VS_ITER);
    }

    /* Vp/Vs-Ratio */
    if (get_float_from_objectlist
        ("VP_VS_RATIO", number_readobjects, &(vinv->VP_VS_RATIO), varname_list, value_list, used_list)) {
        vinv->VP_VS_RATIO = 0.0;
        log_warn("Variable VP_VS_RATIO is set to default value %4.2f which means that it is disregarded.\n",
                 vinv->VP_VS_RATIO);
    }

    /* Limited update of model parameters in reference to the starting model */
    if (get_int_from_objectlist("S", number_readobjects, &(vinv->S), varname_list, value_list, used_list)) {
        vinv->S = 0;
        log_warn("Variable S is set to default value %d.\n", vinv->S);
    } else {
        if (vinv->S == 1) {
            if (get_float_from_objectlist
                ("S_VS", number_readobjects, &(vinv->S_VS), varname_list, value_list, used_list))
                log_fatal("Variable S_VS could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("S_VP", number_readobjects, &(vinv->S_VP), varname_list, value_list, used_list))
                log_fatal("Variable S_VP could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("S_RHO", number_readobjects, &(vinv->S_RHO), varname_list, value_list, used_list))
                log_fatal("Variable S_RHO could not be retrieved from the json input file!");
        }
    }

    /* Upper and lower limits for model parameters */
    if (get_float_from_objectlist
        ("VPUPPERLIM", number_readobjects, &(vinv->VPUPPERLIM), varname_list, value_list, used_list)) {
        vinv->VPUPPERLIM = 25000.0;
        log_warn("Variable VPUPPERLIM is set to default value %f.\n", vinv->VPUPPERLIM);
    }
    if (get_float_from_objectlist
        ("VPLOWERLIM", number_readobjects, &(vinv->VPLOWERLIM), varname_list, value_list, used_list)) {
        vinv->VPLOWERLIM = 0.0;
        log_warn("Variable VPLOWERLIM is set to default value %f.\n", vinv->VPLOWERLIM);
    }
    if (get_float_from_objectlist
        ("VSUPPERLIM", number_readobjects, &(vinv->VSUPPERLIM), varname_list, value_list, used_list)) {
        vinv->VSUPPERLIM = 25000.0;
        log_warn("Variable VSUPPERLIM is set to default value %f.\n", vinv->VSUPPERLIM);
    }
    if (get_float_from_objectlist
        ("VSLOWERLIM", number_readobjects, &(vinv->VSLOWERLIM), varname_list, value_list, used_list)) {
        vinv->VSLOWERLIM = 0.0;
        log_warn("Variable VSLOWERLIM is set to default value %f.\n", vinv->VSLOWERLIM);
    }
    if (get_float_from_objectlist
        ("RHOUPPERLIM", number_readobjects, &(vinv->RHOUPPERLIM), varname_list, value_list, used_list)) {
        vinv->RHOUPPERLIM = 25000.0;
        log_warn("Variable RHOUPPERLIM is set to default value %f.\n", vinv->RHOUPPERLIM);
    }
    if (get_float_from_objectlist
        ("RHOLOWERLIM", number_readobjects, &(vinv->RHOLOWERLIM), varname_list, value_list, used_list)) {
        vinv->RHOLOWERLIM = 0.0;
        log_warn("Variable RHOLOWERLIM is set to default value %f.\n", vinv->RHOLOWERLIM);
    }

    /* Hessian and Gradient-Method */
    if (get_int_from_objectlist
        ("GRAD_METHOD", number_readobjects, &(vinv->GRAD_METHOD), varname_list, value_list, used_list)) {
        log_fatal("Variable GRAD_METHOD could not be retrieved from the json input file!");
    } else {
        if (vinv->GRAD_METHOD == 1) {
            vinv->WOLFE_CONDITION = 0;
        }
        if (vinv->GRAD_METHOD == 2) {
            if (get_int_from_objectlist
                ("LBFGS_STEP_LENGTH", number_readobjects, &(vinv->LBFGS_STEP_LENGTH), varname_list, value_list,
                 used_list)) {
                vinv->LBFGS_STEP_LENGTH = 1;
                log_warn("Variable LBFGS_STEP_LENGTH is set to default value %d.\n", vinv->LBFGS_STEP_LENGTH);
            }
            if (get_int_from_objectlist
                ("N_LBFGS", number_readobjects, &(vinv->N_LBFGS), varname_list, value_list, used_list)) {
                vinv->N_LBFGS = 5;
                log_warn("Variable N_LBFGS is set to default value %d.\n", vinv->N_LBFGS);
            }
            if (get_float_from_objectlist
                ("LBFGS_SCALE_GRADIENTS", number_readobjects, &(vinv->LBFGS_SCALE_GRADIENTS), varname_list, value_list,
                 used_list)) {
                vinv->LBFGS_SCALE_GRADIENTS = 1;
            }
            if (get_int_from_objectlist
                ("WOLFE_CONDITION", number_readobjects, &(vinv->WOLFE_CONDITION), varname_list, value_list,
                 used_list)) {
                vinv->WOLFE_CONDITION = 1;
                log_warn("Variable WOLFE_CONDITION is set to default value %d.\n", vinv->WOLFE_CONDITION);
            } else {
                if (get_int_from_objectlist
                    ("WOLFE_NUM_TEST", number_readobjects, &(vinv->WOLFE_NUM_TEST), varname_list, value_list,
                     used_list)) {
                    vinv->WOLFE_NUM_TEST = 10;
                    log_warn("Variable WOLFE_NUM_TEST is set to default value %d.\n", vinv->WOLFE_NUM_TEST);
                }
                if (get_int_from_objectlist
                    ("WOLFE_TRY_OLD_STEPLENGTH", number_readobjects, &(vinv->WOLFE_TRY_OLD_STEPLENGTH), varname_list,
                     value_list, used_list)) {
                    vinv->WOLFE_TRY_OLD_STEPLENGTH = 0;
                    log_warn("Variable WOLFE_TRY_OLD_STEPLENGTH is set to default value %d.\n",
                             vinv->WOLFE_TRY_OLD_STEPLENGTH);
                }
                if (get_float_from_objectlist
                    ("WOLFE_C1_SL", number_readobjects, &(vinv->WOLFE_C1_SL), varname_list, value_list, used_list)) {
                    vinv->WOLFE_C1_SL = 1e-4;
                    log_warn("Variable WOLFE_C1_SL is set to default value %f.\n", vinv->WOLFE_C1_SL);
                }
                if (get_float_from_objectlist
                    ("WOLFE_C2_SL", number_readobjects, &(vinv->WOLFE_C2_SL), varname_list, value_list, used_list)) {
                    vinv->WOLFE_C2_SL = 0.9;
                    log_warn("Variable WOLFE_C2_SL is set to default value %f.\n", vinv->WOLFE_C2_SL);
                }
            }
        }
    }

    /* Definition of smoothing the models vp and vs */
    if (get_int_from_objectlist
        ("MODEL_FILTER", number_readobjects, &(vinv->MODEL_FILTER), varname_list, value_list, used_list)) {
        vinv->MODEL_FILTER = 0;
        log_warn("Variable MODEL_FILTER is set to default value %d.\n", vinv->MODEL_FILTER);
    } else {
        if (vinv->MODEL_FILTER == 1) {
            if (get_int_from_objectlist
                ("FILT_SIZE", number_readobjects, &(vinv->FILT_SIZE), varname_list, value_list, used_list))
                log_fatal("Variable FILT_SIZE could not be retrieved from the json input file!");
        }
    }

    if (vinv->INV_STF) {
        /* Trace killing STF */
        if (get_int_from_objectlist
            ("TRKILL_STF", number_readobjects, &(vinv->TRKILL_STF), varname_list, value_list, used_list)) {
            vinv->TRKILL_STF = 0;
            log_warn("Variable TRKILL_STF is set to default value %d.\n", vinv->TRKILL_STF);
        } else {
            if (get_int_from_objectlist
                ("STF_FULL", number_readobjects, &(vinv->STF_FULL), varname_list, value_list, used_list)) {
                vinv->STF_FULL = 0;
                log_warn("Variable STF_FULL is set to default value %d.\n", vinv->STF_FULL);
            }
            if (vinv->TRKILL_STF == 1) {
                if (get_int_from_objectlist
                    ("TRKILL_STF_OFFSET", number_readobjects, &(vinv->TRKILL_STF_OFFSET), varname_list, value_list,
                     used_list)) {
                    vinv->TRKILL_STF_OFFSET = 0;
                }
                if (vinv->TRKILL_STF_OFFSET == 0) { /* Only TRKILL File */
                    if (get_string_from_objectlist
                        ("TRKILL_FILE_STF", number_readobjects, (vinv->TRKILL_FILE_STF), varname_list, value_list,
                         used_list))
                        log_fatal("Variable TRKILL_FILE_STF could not be retrieved from the json input file!");
                }
                if (vinv->TRKILL_STF_OFFSET == 1) { /* Only Offset based TRKill */
                    if (get_int_from_objectlist
                        ("TRKILL_STF_OFFSET_INVERT", number_readobjects, &(vinv->TRKILL_STF_OFFSET_INVERT),
                         varname_list, value_list, used_list)) {
                        vinv->TRKILL_STF_OFFSET_INVERT = 0;
                    }
                    if (get_float_from_objectlist
                        ("TRKILL_STF_OFFSET_LOWER", number_readobjects, &(vinv->TRKILL_STF_OFFSET_LOWER), varname_list,
                         value_list, used_list)) {
                        vinv->TRKILL_STF_OFFSET_LOWER = 0.0;
                    }
                    if (get_float_from_objectlist
                        ("TRKILL_STF_OFFSET_UPPER", number_readobjects, &(vinv->TRKILL_STF_OFFSET_UPPER), varname_list,
                         value_list, used_list)) {
                        log_fatal("Variable TRKILL_STF_OFFSET_UPPER could not be retrieved from the json input file!");
                    }
                }
                if (vinv->TRKILL_STF_OFFSET == 2) { /* Both Offset based TRKill & File */
                    if (get_int_from_objectlist
                        ("TRKILL_STF_OFFSET_INVERT", number_readobjects, &(vinv->TRKILL_STF_OFFSET_INVERT),
                         varname_list, value_list, used_list)) {
                        vinv->TRKILL_STF_OFFSET_INVERT = 0;
                    } else {
                        if (vinv->TRKILL_STF_OFFSET_INVERT == 1) {
                            log_fatal("Variable TRKILL_STF_OFFSET_INVERT==1 and TRKILL_STF_OFFSET==2 not possible!");
                        }
                    }
                    if (get_int_from_objectlist
                        ("TRKILL_STF_OFFSET_INVERT", number_readobjects, &(vinv->TRKILL_STF_OFFSET_INVERT),
                         varname_list, value_list, used_list)) {
                        vinv->TRKILL_STF_OFFSET_INVERT = 0;
                    }

                    if (get_float_from_objectlist
                        ("TRKILL_STF_OFFSET_LOWER", number_readobjects, &(vinv->TRKILL_STF_OFFSET_LOWER), varname_list,
                         value_list, used_list)) {
                        vinv->TRKILL_STF_OFFSET_LOWER = 0.0;
                    }
                    if (get_float_from_objectlist
                        ("TRKILL_STF_OFFSET_UPPER", number_readobjects, &(vinv->TRKILL_STF_OFFSET_UPPER), varname_list,
                         value_list, used_list)) {
                        log_fatal("Variable TRKILL_STF_OFFSET_UPPER could not be retrieved from the json input file!");
                    }
                    if (get_string_from_objectlist
                        ("TRKILL_FILE_STF", number_readobjects, (vinv->TRKILL_FILE_STF), varname_list, value_list,
                         used_list))
                        log_fatal("Variable TRKILL_FILE_STF could not be retrieved from the json input file!");
                }
            }
        }
        /* Taper STF */
        if (get_int_from_objectlist
            ("TAPER_STF", number_readobjects, &(vinv->TAPER_STF), varname_list, value_list, used_list)) {
            vinv->TAPER_STF = 0;
            log_warn("Variable TAPER_STF is set to default value %d.\n", vinv->TAPER_STF);
        }
    }

    /* Frequency filtering during inversion */
    if (get_int_from_objectlist
        ("TIME_FILT", number_readobjects, &(vinv->TIME_FILT), varname_list, value_list, used_list)) {
        vinv->TIME_FILT = 0;
        log_warn("Variable TIME_FILT is set to default value %d.\n", vinv->TIME_FILT);
    } else {
        if (get_int_from_objectlist
            ("WRITE_FILTERED_DATA", number_readobjects, &(vinv->WRITE_FILTERED_DATA), varname_list, value_list,
             used_list)) {
            vinv->WRITE_FILTERED_DATA = 0;
        }
        if (get_float_from_objectlist
            ("F_HIGH_PASS", number_readobjects, &(vinv->F_HIGH_PASS), varname_list, value_list, used_list)) {
            vinv->F_HIGH_PASS = 0.0;
            log_warn("Variable F_HIGH_PASS is set to default value %f.\n", vinv->F_HIGH_PASS);
        }
        if (vinv->TIME_FILT == 1) {
            if (get_float_from_objectlist
                ("F_LOW_PASS_START", number_readobjects, &(vinv->F_LOW_PASS_START), varname_list, value_list,used_list)) {
                log_fatal("Variable F_LOW_PASS_START could not be retrieved from the json input file!");
            }
            if (get_float_from_objectlist
                ("F_LOW_PASS_END", number_readobjects, &(vinv->F_LOW_PASS_END), varname_list, value_list, used_list)) {
                log_fatal("Variable F_LOW_PASS_END could not be retrieved from the json input file!");
            }
            if (get_float_from_objectlist
                ("F_LOW_PASS_INCR", number_readobjects, &(vinv->F_LOW_PASS_INCR), varname_list, value_list,
                 used_list)) {
                log_fatal("Variable F_LOW_PASS_INCR could not be retrieved from the json input file!");
            }
            if (get_int_from_objectlist
                ("ORDER", number_readobjects, &(vinv->ORDER), varname_list, value_list, used_list)) {
                log_fatal("Variable ORDER could not be retrieved from the json input file!");
            }
        }
        if (vinv->TIME_FILT == 2) {
            if (get_float_from_objectlist
                ("F_HIGH_PASS", number_readobjects, &(vinv->F_HIGH_PASS), varname_list, value_list, used_list)) {
                vinv->F_HIGH_PASS = 0.0;
                log_warn("Variable F_HIGH_PASS is set to default value %f.\n", vinv->F_HIGH_PASS);
            }
            if (get_string_from_objectlist
                ("FREQ_FILE", number_readobjects, (vinv->FREQ_FILE), varname_list, value_list, used_list))
                log_fatal("Variable FREQ_FILE could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("ORDER", number_readobjects, &(vinv->ORDER), varname_list, value_list, used_list))
                log_fatal("Variable ORDER could not be retrieved from the json input file!");
        }
    }

    /* Gradient calculation */
    if (get_int_from_objectlist("LNORM", number_readobjects, &(vinv->LNORM), varname_list, value_list, used_list))
        log_fatal("Variable LNORM could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("NORMALIZE", number_readobjects, &(vinv->NORMALIZE), varname_list, value_list, used_list))
        log_fatal("Variable NORMALIZE could not be retrieved from the json input file!");
    if (get_int_from_objectlist("DTINV", number_readobjects, &(vinv->DTINV), varname_list, value_list, used_list)) {
        log_fatal("Variable DTINV could not be retrieved from the json input file!");
    } else {
        if (vinv->DTINV <= 0) log_fatal("Variable DTINV must be a positive integer!");
        vinv->NTDTINV = ceil((float)gv->NT / (float)vinv->DTINV);
    }
    if (vinv->LNORM == 8) {
        if (get_float_from_objectlist
            ("WATERLEVEL_LNORM8", number_readobjects, &(vinv->WATERLEVEL_LNORM8), varname_list, value_list, used_list))
            log_fatal("Variable WATERLEVEL_LNORM8 could not be retrieved from the json input file!");
    }

    /* Step length estimation */
    if (get_float_from_objectlist
        ("EPS_SCALE", number_readobjects, &(vinv->EPS_SCALE), varname_list, value_list, used_list))
        log_fatal("Variable EPS_SCALE could not be retrieved from the json input file!");
    if (get_int_from_objectlist("STEPMAX", number_readobjects, &(vinv->STEPMAX), varname_list, value_list, used_list))
        log_fatal("Variable STEPMAX could not be retrieved from the json input file!");
    if (get_float_from_objectlist
        ("SCALEFAC", number_readobjects, &(vinv->SCALEFAC), varname_list, value_list, used_list))
        log_fatal("Variable SCALEFAC could not be retrieved from the json input file!");

    /* Termination of the program */
    if (get_float_from_objectlist("PRO", number_readobjects, &(vinv->PRO), varname_list, value_list, used_list))
        log_fatal("Variable PRO could not be retrieved from the json input file!");
    if (get_int_from_objectlist("MIN_ITER", number_readobjects, &(vinv->MIN_ITER), varname_list, value_list, used_list)) {
        vinv->MIN_ITER = 0;
        log_warn("Variable MIN_ITER is set to default value %d.\n", vinv->MIN_ITER);
    }

    /* Trace killing */
    if (get_int_from_objectlist("TRKILL", number_readobjects, &(vinv->TRKILL), varname_list, value_list, used_list)) {
        vinv->TRKILL = 0;
        log_warn("Variable TRKILL is set to default value %d.\n", vinv->TRKILL);
    } else {
        if (vinv->TRKILL == 1) {
            if (get_int_from_objectlist
                ("TRKILL_OFFSET", number_readobjects, &(vinv->TRKILL_OFFSET), varname_list, value_list, used_list)) {
                vinv->TRKILL_OFFSET = 0;
            }
            if (vinv->TRKILL_OFFSET == 0) { /* Only TRKILL File */
                if (get_string_from_objectlist
                    ("TRKILL_FILE", number_readobjects, (vinv->TRKILL_FILE), varname_list, value_list, used_list))
                    log_fatal("Variable TRKILL_FILE could not be retrieved from the json input file!");
            } else if (vinv->TRKILL_OFFSET == 1) {  /* Only Offset based TRKill */
                if (get_float_from_objectlist
                    ("TRKILL_OFFSET_LOWER", number_readobjects, &(vinv->TRKILL_OFFSET_LOWER), varname_list, value_list,
                     used_list)) {
                    vinv->TRKILL_OFFSET_LOWER = 0.0;
                }
                if (get_float_from_objectlist
                    ("TRKILL_OFFSET_UPPER", number_readobjects, &(vinv->TRKILL_OFFSET_UPPER), varname_list, value_list,
                     used_list)) {
                    log_fatal("Variable TRKILL_OFFSET_UPPER could not be retrieved from the json input file!");
                }
            } else if (vinv->TRKILL_OFFSET == 2) {  /* Both Offset based TRKill & File */
                if (get_float_from_objectlist
                    ("TRKILL_OFFSET_LOWER", number_readobjects, &(vinv->TRKILL_OFFSET_LOWER), varname_list, value_list,
                     used_list)) {
                    vinv->TRKILL_OFFSET_LOWER = 0.0;
                }
                if (get_float_from_objectlist
                    ("TRKILL_OFFSET_UPPER", number_readobjects, &(vinv->TRKILL_OFFSET_UPPER), varname_list, value_list,
                     used_list)) {
                    log_fatal("Variable TRKILL_OFFSET_UPPER could not be retrieved from the json input file!");
                }
                if (get_string_from_objectlist
                    ("TRKILL_FILE", number_readobjects, (vinv->TRKILL_FILE), varname_list, value_list, used_list))
                    log_fatal("Variable TRKILL_FILE could not be retrieved from the json input file!");
            } else {
                log_fatal("Variable TRKILL_OFFSET can only have the values 0, 1 or 2!");
            }
        }

        /* Time windowing and VPPML */
        if (get_int_from_objectlist
            ("TIMEWIN", number_readobjects, &(vinv->TIMEWIN), varname_list, value_list, used_list)) {
            vinv->TIMEWIN = 0;
            log_warn("Variable TIMEWIN is set to default value %d.\n", vinv->TIMEWIN);
        } else {
            if (vinv->TIMEWIN == 1) {
                if (get_int_from_objectlist
                    ("TW_IND", number_readobjects, &(vinv->TW_IND), varname_list, value_list, used_list)) {
                    vinv->TW_IND = 0;
                    log_warn("Variable TW_IND is set to default value %d.\n", vinv->TW_IND);
                } else {
                    if (vinv->TW_IND > 2)
                        log_fatal("Only TW_IND=1 (one time window) and TW_IND=2 (two time windows) possible");
                }
                if (get_string_from_objectlist
                    ("PICKS_FILE", number_readobjects, (vinv->PICKS_FILE), varname_list, value_list, used_list))
                    log_fatal("Variable PICKS_FILE could not be retrieved from the json input file!");
                if (get_float_from_objectlist
                    ("TWLENGTH_PLUS", number_readobjects, &(vinv->TWLENGTH_PLUS), varname_list, value_list, used_list))
                    log_fatal("Variable TWLENGTH_PLUS could not be retrieved from the json input file!");
                if (get_float_from_objectlist
                    ("TWLENGTH_MINUS", number_readobjects, &(vinv->TWLENGTH_MINUS), varname_list, value_list,
                     used_list))
                    log_fatal("Variable TWLENGTH_MINUS could not be retrieved from the json input file!");
                if (get_float_from_objectlist
                    ("GAMMA", number_readobjects, &(vinv->GAMMA), varname_list, value_list, used_list))
                    log_fatal("Variable GAMMA could not be retrieved from the json input file!");
            }
        }
    }

    /********************************************/
    /* Check files and directories if necessary */

    /********************************************/

    /* workflow file */
    if (vinv->USE_WORKFLOW != 0) {
        if (access((vinv->FILE_WORKFLOW), 0) != 0) {
            log_error("Workflow file %s does not exists.\n", (vinv->FILE_WORKFLOW));
            fserr = 1;
        } else if (access((vinv->FILE_WORKFLOW), 4) != 0) {
            log_error("Workflow file %s cannot be read, no read permission.\n", (vinv->FILE_WORKFLOW));
            fserr = 1;
        }
    }

    /* frequency file */
    if (vinv->TIME_FILT == 2) {
        if (access((vinv->FREQ_FILE), 0) != 0) {
            log_error("Frequency file %s does not exists.\n", (vinv->FREQ_FILE));
            fserr = 1;
        } else if (access((vinv->FREQ_FILE), 4) != 0) {
            log_error("Frequency file %s cannot be read, no read permission.\n", (vinv->FREQ_FILE));
            fserr = 1;
        }
    }

    /* picks file */
    if (vinv->TIMEWIN == 1) {
        if (access((vinv->PICKS_FILE), 0) != 0) {
            log_error("Picks file %s does not exists.\n", (vinv->PICKS_FILE));
            fserr = 1;
        } else if (access((vinv->PICKS_FILE), 4) != 0) {
            log_error("Picks file %s cannot be read, no read permission.\n", (vinv->PICKS_FILE));
            fserr = 1;
        }
    }

    /* trace kill file */
    if (vinv->TRKILL == 1) {
        if (access((vinv->TRKILL_FILE), 0) != 0) {
            log_error("Trace kill file %s does not exists.\n", (vinv->TRKILL_FILE));
            fserr = 1;
        } else if (access((vinv->TRKILL_FILE), 4) != 0) {
            log_error("Trace kill file %s cannot be read, no read permission.\n", (vinv->TRKILL_FILE));
            fserr = 1;
        }
    }

    return fserr;

}
