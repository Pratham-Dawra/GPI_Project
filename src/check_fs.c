
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

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>

#include "fd.h"
#include "globvar_struct.h"
#include "logging.h"

void check_fs(GlobVar * gv)
{
    int fserr = 0;

  /********************************************/
    /* Check output directories as required     */

  /********************************************/

    if (gv->LOG > 0) {
        /* check log file directory */
        char *dirc = strdup(gv->LOG_FILE);
        if (!dirc)
            log_fatal("Could not copy string in check_fs.c - memory issue encountered.\n");
        char *dname = dirname(dirc);
        if (access(dname, W_OK) != 0) {
            log_error("Cannot write to log directory %s.\n", dname);
            fserr += 1;
        } else {
            log_infoc(0, "Filesystem check: log file directory %s is writable.\n", dname);
        }
        free(dirc);
    }

    if (gv->SEISMO > 0) {
        /* check seismogram directory */
        char *dirc = strdup(gv->SEIS_FILE);
        if (!dirc)
            log_fatal("Could not copy string in check_fs.c - memory issue encountered.\n");
        char *dname = dirname(dirc);
        if (access(dname, W_OK) != 0) {
            log_error("Cannot write to seismogram directory %s.\n", dname);
            fserr += 1;
        } else {
            log_infoc(0, "Filesystem check: seismogram directory %s is writable.\n", dname);
        }
        free(dirc);
    }

    if (gv->SNAP > 0) {
        /* check snapshot directory */
        char *dirc = strdup(gv->SNAP_FILE);
        if (!dirc)
            log_fatal("Could not copy string in check_fs.c - memory issue encountered.\n");
        char *dname = dirname(dirc);
        if (access(dname, W_OK) != 0) {
            log_error("Cannot write to snapshot directory %s.\n", dname);
            fserr += 1;
        } else {
            log_infoc(0, "Filesystem check: snapshot directory %s is writable.\n", dname);
        }
        free(dirc);
    }

    if (gv->WRITE_MODELFILES > 0) {
        /* check model file directory */
        char *dirc = strdup(gv->MFILE);
        if (!dirc)
            log_fatal("Could not copy string in check_fs.c - memory issue encountered.\n");
        char *dname = dirname(dirc);
        if (access(dname, W_OK) != 0) {
            log_error("Cannot write to model file directory %s.\n", dname);
            fserr += 1;
        } else {
            log_infoc(0, "Filesystem check: model file directory %s is writable.\n", dname);
        }
        free(dirc);
    }

    if (1 == gv->SIGOUT) {
        /* check signal output directory */
        char *dirc = strdup(gv->SIGOUT_FILE);
        if (!dirc)
            log_fatal("Could not copy string in check_fs.c - memory issue encountered.\n");
        char *dname = dirname(dirc);
        if (access(dname, W_OK) != 0) {
            log_error("Cannot write to source output directory %s.\n", dname);
            fserr += 1;
        } else {
            log_infoc(0, "Filesystem check: source output directory %s is writable.\n", dname);
        }
        free(dirc);
    }

  /********************************************/
    /* ERROR                                    */

  /********************************************/
    if (fserr>0) {
        log_fatal("%d error(s) encountered while checking filesystem for correct permissions.\n",fserr);
    }

    return;
}
