
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
 *   program SOFI2D, reading input-parameters from input-file or stdin
 *   that are formatted according to the json standard
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include "fd.h"
#include "logging.h"
#include "enums.h"

static const char *weq_descr[NWEQ] = { "AC_ISO", "AC_VTI", "AC_TTI",
    "EL_ISO", "VEL_ISO", "EL_VTI",
    "VEL_VTI", "EL_TTI", "VEL_TTI",
    "VAC_ISO", "VAC_VTI", "VAC_TTI"
};

static const char *weq_verbose[NWEQ] = { "acoustic wave equation",
    "acoustic VTI wave equation",
    "acoustic TTI wave equation",
    "elastic wave equation",
    "viscoelastic wave equation",
    "elastic VTI wave equation",
    "viscoelastic VTI wave equation",
    "elastic TTI wave equation",
    "viscoelastic TTI wave equation",
    "viscoacoustic wave equation",
    "viscoacoustic VTI wave equation",
    "viscoacoustic TTI wave equation"
};

const char *get_weq_verbose(WEQTYPE wt)
{
    return weq_verbose[wt];
}

static void parse_weqtype(const char *weqt, GlobVar * gv)
{
    bool b_found = false;

    for (int i = 0; i < NWEQ; ++i) {
        if (!strncasecmp(weqt, weq_descr[i], STRING_SIZE - 1)) {
            gv->WEQ = (WEQTYPE) i;
            b_found = true;
            break;
        }
    }

    if (!b_found) {
        log_error("Unknown value '%s' for parameter WEQ.\n", weqt);
        log_fatal("Unknown wave equation (WEQ) specified in parameter file.\n");
    }
}

void read_par_json(const char *fileinp, GlobVar * gv)
{
    /* definition of local variables */
    int number_readobjects = 0, fserr = 0, i;
    int *used_list = NULL, not_used = 0;
    bool header_printed = false;

    //allocate first object in list
    char **varname_list = malloc(STRING_SIZE * sizeof(char *));
    char **value_list = malloc(STRING_SIZE * sizeof(char *));

    //read in objects from file
    number_readobjects = read_objects_from_intputfile(fileinp, varname_list, value_list);
    log_info("Read %d parameters from file %s.\n", number_readobjects, fileinp);

    used_list = (int *)calloc(sizeof(int), number_readobjects);

    if (get_string_from_objectlist
        ("LOG_VERBOSITY", number_readobjects, (gv->LOG_VERBOSITY), varname_list, value_list, used_list)) {
        strcpy(&(gv->LOG_VERBOSITY[0]), "INFO");
        log_set_level(LOG_INFO);
        log_warn("Parameter LOG_VERBOSITY not given, set to 'INFO' (default)!\n");
    }
    log_set_level_from_string(&(gv->LOG_VERBOSITY[0]));

    //print objects to screen
    print_objectlist_screen(number_readobjects, varname_list, value_list);

    // extract variables from object list

    /* =================================
     * section general grid and discretization parameters
     * ================================= */

    if (get_int_from_objectlist("NPROCX", number_readobjects, &(gv->NPROCX), varname_list, value_list, used_list))
        log_fatal("Variable NPROCX could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NPROCY", number_readobjects, &(gv->NPROCY), varname_list, value_list, used_list))
        log_fatal("Variable NPROCY could not be retrieved from the json input file!");
    gv->NPROC = (gv->NPROCX) * (gv->NPROCY);
    /*    if (get_int_from_objectlist("RSG",number_readobjects,&(gv->RSG),varname_list, value_list, used_list))
     * log_fatal("Variable RSG could not be retrieved from the json input file!"); */
    if (get_int_from_objectlist("FDORDER", number_readobjects, &(gv->FDORDER), varname_list, value_list, used_list))
        log_fatal("Variable FDORDER could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("FDORDER_TIME", number_readobjects, &(gv->FDORDER_TIME), varname_list, value_list, used_list)) {
        (gv->FDORDER_TIME) = 2;
    } else {
        if ((gv->FDORDER_TIME) != 2 && (gv->FDORDER_TIME) != 4) {
            log_fatal("Only FDORDER_TIME 2 or 4 are supported!");
        }
    }
    if (get_int_from_objectlist
        ("MAXRELERROR", number_readobjects, &(gv->MAXRELERROR), varname_list, value_list, used_list))
        log_fatal("Variable MAXRELERROR could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NX", number_readobjects, &(gv->NXG), varname_list, value_list, used_list))
        log_fatal("Variable NX could not be retrieved from the json input file!");
    if (get_int_from_objectlist("NY", number_readobjects, &(gv->NYG), varname_list, value_list, used_list))
        log_fatal("Variable NY could not be retrieved from the json input file!");
    if (get_float_from_objectlist("DH", number_readobjects, &(gv->DH), varname_list, value_list, used_list))
        log_fatal("Variable DH could not be retrieved from the json input file!");
    if (get_float_from_objectlist("TIME", number_readobjects, &(gv->TIME), varname_list, value_list, used_list))
        log_fatal("Variable TIME could not be retrieved from the json input file!");
    if (get_float_from_objectlist("DT", number_readobjects, &(gv->DT), varname_list, value_list, used_list))
        log_fatal("Variable DT could not be retrieved from the json input file!");

    /* =================================
     * wave equation
     * ================================= */

    char weqtype[STRING_SIZE];

    if (get_string_from_objectlist("WEQ", number_readobjects, weqtype, varname_list, value_list, used_list)) {
        gv->WEQ = EL_ISO;
    } else {
        parse_weqtype(weqtype, gv);
    }

    /* =================================
     * section source parameters
     * ================================= */

    if (get_int_from_objectlist
        ("SOURCE_TYPE", number_readobjects, &(gv->SOURCE_TYPE), varname_list, value_list, used_list))
        log_fatal("Variable SOURCE_TYPE could not be retrieved from the json input file!");
    if (get_int_from_objectlist("SIGOUT", number_readobjects, &(gv->SIGOUT), varname_list, value_list, used_list)) {
        log_info
            ("For your information: variables SIGOUT, SIGOUT_FILE & SIGOUT_FORMAT available to output source signal to file.\n");
    } else {
        if ((gv->SIGOUT) == 1) {
            if (get_string_from_objectlist
                ("SIGOUT_FILE", number_readobjects, gv->SIGOUT_FILE, varname_list, value_list, used_list)) {
                log_fatal("Variable SIGOUT_FILE could not be retrieved from the json input file!");
            }
            if (get_int_from_objectlist
                ("SIGOUT_FORMAT", number_readobjects, &(gv->SIGOUT_FORMAT), varname_list, value_list, used_list)) {
                log_fatal("Variable SIGOUT_FORMAT could not be retrieved from the json input file!");
            }
        }
    }
    if (get_int_from_objectlist
        ("SOURCE_SHAPE", number_readobjects, &(gv->SOURCE_SHAPE), varname_list, value_list, used_list))
        log_fatal("Variable SOURCE_SHAPE could not be retrieved from the json input file!");
    else {
        if ((gv->SOURCE_SHAPE) == 3) {
            if (get_string_from_objectlist
                ("SIGNAL_FILE", number_readobjects, (gv->SIGNAL_FILE), varname_list, value_list, used_list))
                log_fatal("Variable SIGNAL_FILE could not be retrieved from the json input file!");
        }
    }
    if (get_int_from_objectlist("SRCREC", number_readobjects, &(gv->SRCREC), varname_list, value_list, used_list))
        log_fatal("Variable SRCREC could not be retrieved from the json input file!");
    else {
        if ((gv->SRCREC) == 1) {
            if (get_string_from_objectlist
                ("SOURCE_FILE", number_readobjects, (gv->SOURCE_FILE), varname_list, value_list, used_list))
                log_fatal("Variable SOURCE_FILE could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("RUN_MULTIPLE_SHOTS", number_readobjects, &(gv->RUN_MULTIPLE_SHOTS), varname_list, value_list,
                 used_list))
                log_fatal("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
        }
        if ((gv->SRCREC) == 2) {
            if (get_float_from_objectlist
                ("PLANE_WAVE_DEPTH", number_readobjects, &(gv->PLANE_WAVE_DEPTH), varname_list, value_list, used_list))
                log_fatal("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
            else {
                if ((gv->PLANE_WAVE_DEPTH) > 0.0) {
                    if (get_float_from_objectlist
                        ("PLANE_WAVE_ANGLE", number_readobjects, &(gv->PLANE_WAVE_ANGLE), varname_list, value_list,
                         used_list))
                        log_fatal("Variable PLANE_WAVE_ANGLE could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("TS", number_readobjects, &(gv->TS), varname_list, value_list, used_list))
                        log_fatal("Variable TS could not be retrieved from the json input file!");
                }
            }
        }
    }

    /* =================================
     * section boundary parameters
     * ================================= */

    if (get_int_from_objectlist("FREE_SURF", number_readobjects, &(gv->FREE_SURF), varname_list, value_list, used_list))
        log_fatal("Variable FREE_SURF could not be retrieved from the json input file!");
    if (get_int_from_objectlist("BOUNDARY", number_readobjects, &(gv->BOUNDARY), varname_list, value_list, used_list))
        log_fatal("Variable BOUNDARY could not be retrieved from the json input file!");
    if (get_int_from_objectlist("FW", number_readobjects, &(gv->FW), varname_list, value_list, used_list))
        log_fatal("Variable FW could not be retrieved from the json input file!");
    if (get_int_from_objectlist("ABS_TYPE", number_readobjects, &(gv->ABS_TYPE), varname_list, value_list, used_list))
        log_fatal("Variable ABS_TYPE could not be retrieved from the json input file!");

    if ((gv->ABS_TYPE) == 1) {
        if (get_float_from_objectlist("NPOWER", number_readobjects, &(gv->NPOWER), varname_list, value_list, used_list))
            log_fatal("Variable NPOWER could not be retrieved from the json input file!");
        if (get_float_from_objectlist
            ("K_MAX_CPML", number_readobjects, &(gv->K_MAX_CPML), varname_list, value_list, used_list))
            log_fatal("Variable K_MAX_CPML could not be retrieved from the json input file!");
        if (get_float_from_objectlist("FPML", number_readobjects, &(gv->FPML), varname_list, value_list, used_list))
            log_fatal("Variable FPML could not be retrieved from the json input file!");
        if (get_float_from_objectlist("VPPML", number_readobjects, &(gv->VPPML), varname_list, value_list, used_list))
            log_fatal("Variable VPPML could not be retrieved from the json input file!");
    }
    if ((gv->ABS_TYPE) == 2) {
        if (get_float_from_objectlist
            ("DAMPING", number_readobjects, &(gv->DAMPING), varname_list, value_list, used_list))
            log_fatal("Variable DAMPING could not be retrieved from the json input file!");
    }

    /* =================================
     * section snapshot parameters
     * ================================= */

    if (get_int_from_objectlist("SNAP", number_readobjects, &(gv->SNAP), varname_list, value_list, used_list))
        log_fatal("Variable SNAP could not be retrieved from the json input file!");
    else {
        if ((gv->SNAP) > 0) {
            if (get_int_from_objectlist
                ("SNAP_FORMAT", number_readobjects, &(gv->SNAP_FORMAT), varname_list, value_list, used_list))
                log_fatal("Variable SNAP_FORMAT could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("TSNAP1", number_readobjects, &(gv->TSNAP1), varname_list, value_list, used_list))
                log_fatal("Variable TSNAP1 could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("TSNAP2", number_readobjects, &(gv->TSNAP2), varname_list, value_list, used_list))
                log_fatal("Variable TSNAP2 could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("TSNAPINC", number_readobjects, &(gv->TSNAPINC), varname_list, value_list, used_list))
                log_fatal("Variable TSNAPINC could not be retrieved from the json input file!");
            if (get_string_from_objectlist
                ("SNAP_FILE", number_readobjects, (gv->SNAP_FILE), varname_list, value_list, used_list))
                log_fatal("Variable SNAP_FILE could not be retrieved from the json input file!");
        }
    }
    /* increments are read in any case, because they will be also used as increment for model output */
    if (get_int_from_objectlist("IDX", number_readobjects, &(gv->IDX), varname_list, value_list, used_list))
        log_fatal("Variable IDX could not be retrieved from the json input file!");
    if (get_int_from_objectlist("IDY", number_readobjects, &(gv->IDY), varname_list, value_list, used_list))
        log_fatal("Variable IDY could not be retrieved from the json input file!");

    /* =================================
     * section seismogramm parameters
     * ================================= */

    if (get_int_from_objectlist("SEISMO", number_readobjects, &(gv->SEISMO), varname_list, value_list, used_list))
        log_fatal("Variable SEISMO could not be retrieved from the json input file!");
    else {
        if ((gv->SEISMO) > 0) {
            if (get_string_from_objectlist
                ("REC_FILE", number_readobjects, (gv->REC_FILE), varname_list, value_list, used_list))
                log_fatal("Variable REC_FILE could not be retrieved from the json input file!");
            if (get_string_from_objectlist
                ("SEIS_FILE", number_readobjects, (gv->SEIS_FILE), varname_list, value_list, used_list))
                log_fatal("Variable SEIS_FILE could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("READREC", number_readobjects, &(gv->READREC), varname_list, value_list, used_list))
                log_fatal("Variable READREC could not be retrieved from the json input file!");
            else {
                if ((gv->READREC) == 0) {
                    if (get_float_from_objectlist
                        ("XREC1", number_readobjects, &(gv->XREC1), varname_list, value_list, used_list))
                        log_fatal("Variable XREC1 could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("XREC2", number_readobjects, &(gv->XREC2), varname_list, value_list, used_list))
                        log_fatal("Variable XREC2T could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("YREC1", number_readobjects, &(gv->YREC1), varname_list, value_list, used_list))
                        log_fatal("Variable YREC1 could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("YREC2", number_readobjects, &(gv->YREC2), varname_list, value_list, used_list))
                        log_fatal("Variable YREC2 could not be retrieved from the json input file!");
                }
            }
            if (get_int_from_objectlist("NDT", number_readobjects, &(gv->NDT), varname_list, value_list, used_list))
                log_fatal("Variable NDT could not be retrieved from the json input file!");
            if (get_int_from_objectlist
                ("SEIS_FORMAT", number_readobjects, &(gv->SEIS_FORMAT), varname_list, value_list, used_list))
                log_fatal("Variable SEIS_FORMAT could not be retrieved from the json input file!");

            if (get_int_from_objectlist
                ("REC_ARRAY", number_readobjects, &(gv->REC_ARRAY), varname_list, value_list, used_list))
                log_fatal("Variable REC_ARRAY could not be retrieved from the json input file!");
            else {
                if ((gv->REC_ARRAY) > 0) {
                    if (get_int_from_objectlist
                        ("DRX", number_readobjects, &(gv->DRX), varname_list, value_list, used_list))
                        log_fatal("Variable DRX could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("REC_ARRAY_DEPTH", number_readobjects, &(gv->REC_ARRAY_DEPTH), varname_list, value_list,
                         used_list))
                        log_fatal("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
                    if (get_float_from_objectlist
                        ("REC_ARRAY_DIST", number_readobjects, &(gv->REC_ARRAY_DIST), varname_list, value_list,
                         used_list))
                        log_fatal("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
                }
            }
            if (get_float_from_objectlist
                ("REFRECX", number_readobjects, &(gv->REFREC[1]), varname_list, value_list, used_list))
                log_fatal("Variable REFRECX could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("REFRECY", number_readobjects, &(gv->REFREC[2]), varname_list, value_list, used_list))
                log_fatal("Variable REFRECY could not be retrieved from the json input file!");
            if (get_float_from_objectlist
                ("NGEOPH", number_readobjects, &(gv->NGEOPH), varname_list, value_list, used_list))
                log_fatal("Variable NGEOPH could not be retrieved from the json input file!");
        }
    }

    /* =========================================
     * section general model and log parameters
     * ========================================= */

    if (get_string_from_objectlist("MFILE", number_readobjects, (gv->MFILE), varname_list, value_list, used_list))
        log_fatal("Variable MFILE could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("WRITE_MODELFILES", number_readobjects, &(gv->WRITE_MODELFILES), varname_list, value_list, used_list)) {
    }
    if (get_int_from_objectlist("LOG", number_readobjects, &(gv->LOG), varname_list, value_list, used_list))
        log_fatal("Variable LOG could not be retrieved from the json input file!");
    if (get_string_from_objectlist("LOG_FILE", number_readobjects, (gv->LOG_FILE), varname_list, value_list, used_list))
        log_fatal("Variable LOG_FILE could not be retrieved from the json input file!");
    if (get_int_from_objectlist("READMOD", number_readobjects, &(gv->READMOD), varname_list, value_list, used_list))
        log_fatal("Variable READMOD could not be retrieved from the json input file!");
    if (get_int_from_objectlist
        ("OUT_TIMESTEP_INFO", number_readobjects, &(gv->OUTNTIMESTEPINFO), varname_list, value_list, used_list)) {
    }

    if (get_float_from_objectlist("TAU", number_readobjects, &(gv->TAU), varname_list, value_list, used_list))
        log_fatal("Variable TAU could not be retrieved from the json input file!");
    if (get_int_from_objectlist("L", number_readobjects, &(gv->L), varname_list, value_list, used_list))
        log_fatal("Variable L could not be retrieved from the json input file!");

    if (((gv->L) == 0) && (((gv->WEQ) == VEL_ISO) || ((gv->WEQ) == VEL_VTI) || ((gv->WEQ) == VEL_TTI) ||
                           ((gv->WEQ) == VAC_ISO) || ((gv->WEQ) == VAC_VTI) || ((gv->WEQ) == VAC_TTI)))
        log_fatal("L>0 required for viscoacoustic/elastic simulation!");

    if (((gv->L) > 0) && (((gv->WEQ) == AC_ISO) || ((gv->WEQ) == AC_VTI) || ((gv->WEQ) == AC_TTI) ||
                          ((gv->WEQ) == EL_ISO) || ((gv->WEQ) == EL_VTI) || ((gv->WEQ) == EL_TTI))) {
        log_warn("L reset to zero for elastic/acoustic simulation.\n");
        (gv->L) = 0;
    }

    /* do NOT remove the FALLTHRU comments below, they are used to tell the compiler
     * that this is an intentional fall through */
    if ((gv->L) > 0) {
        if (get_float_from_objectlist("F_REF", number_readobjects, &(gv->F_REF), varname_list, value_list, used_list))
            log_fatal("Variable F_REF could not be retrieved from the json input file!");
        (gv->FL) = vector(1, (gv->L));
        switch ((gv->L)) {
          case 0:
              break;
          case 5:
              if (get_float_from_objectlist
                  ("FL5", number_readobjects, &(gv->FL[5]), varname_list, value_list, used_list))
                  log_fatal("Variable FL5 could not be retrieved from the json input file!");
              /* FALLTHRU */
          case 4:
              if (get_float_from_objectlist
                  ("FL4", number_readobjects, &(gv->FL[4]), varname_list, value_list, used_list))
                  log_fatal("Variable FL4 could not be retrieved from the json input file!");
              /* FALLTHRU */
          case 3:
              if (get_float_from_objectlist
                  ("FL3", number_readobjects, &(gv->FL[3]), varname_list, value_list, used_list))
                  log_fatal("Variable FL3 could not be retrieved from the json input file!");
              /* FALLTHRU */
          case 2:
              if (get_float_from_objectlist
                  ("FL2", number_readobjects, &(gv->FL[2]), varname_list, value_list, used_list))
                  log_fatal("Variable FL2 could not be retrieved from the json input file!");
              /* FALLTHRU */
          case 1:
              if (get_float_from_objectlist
                  ("FL1", number_readobjects, &(gv->FL[1]), varname_list, value_list, used_list))
                  log_fatal("Variable FL1 could not be retrieved from the json input file!");
              break;
          default:
              log_fatal("More than five relaxation parameters (L>5) not implemented!");
              break;
        }
    } else {
        (gv->FL) = NULL;
    }

    /********************************************/
    /* Check files and directories if necessary */
    /********************************************/

    /* signal file */
    if ((gv->SOURCE_SHAPE) == 3) {
        if (access((gv->SIGNAL_FILE), 0) != 0) {
            log_error("Signal file %s does not exists.\n", (gv->SIGNAL_FILE));
            fserr = 1;
        } else if (access((gv->SIGNAL_FILE), 4) != 0) {
            log_error("Signal file %s cannot be read, no read permission.\n", (gv->SIGNAL_FILE));
            fserr = 1;
        }
    }

    /* source file */
    if ((gv->SRCREC) == 1) {
        if (access((gv->SOURCE_FILE), 0) != 0) {
            log_error("Source file %s does not exist.\n", (gv->SOURCE_FILE));
            fserr = 1;
        } else if (access((gv->SOURCE_FILE), 4) != 0) {
            log_error("Source file %s cannot be read, no read permission.\n", (gv->SOURCE_FILE));
            fserr = 1;
        }
    }

    /* receiver file */
    if ((gv->READREC)) {
        if (access((gv->REC_FILE), 0) != 0) {
            log_error("Receiver file %s does not exist.\n", (gv->REC_FILE));
            fserr = 1;
        } else if (access((gv->REC_FILE), 4) != 0) {
            log_error("Receiver file %s cannot be read, no read permission.\n", (gv->REC_FILE));
            fserr = 1;
        }
    }

    /********************************************/
    /* ERROR                                    */
    /********************************************/
    if (fserr) {
        log_fatal("Error(s) encountered while processing json parameter file.\n");
    }

    for (i = 0; i < number_readobjects; ++i) {
        if (!used_list[i]) {
            if (!header_printed) {
                log_info("The following parameter(s) from the json file were unclaimed in this run:\n");
                header_printed = true;
            }
            log_info("  Parameter '%s' with value '%s'.\n", varname_list[i], value_list[i]);
            ++not_used;
        }
    }
    if (header_printed)
        log_warn("You may want to check %d unclaimed parameter(s).\n", not_used);

    free(used_list);

    return;
}
