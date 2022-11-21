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

void check_fs(FILE *fp, GlobVar *gv) {
  
  char errmsg[80];
  int fserr = 0, MYID;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

  /********************************************/
  /* Check output directories as required     */
  /********************************************/

  if (gv->LOG > 1) {
    /* check log file directory */
    char *dirc = strdup(gv->LOG_FILE);
    if (!dirc) declare_error("Could not copy string in check_fs.c - memory issue encountered.");
    char *dname = dirname(dirc);
    if (access(dname, W_OK) != 0) {
      fprintf(fp, "\n==================================================================\n");
      fprintf(fp, "ERROR: PE=%d, cannot write to log directory %s.\n", MYID, dname);
      fprintf(fp, "\n==================================================================\n");
      fserr = 1;
    } else {
      if (MYID == 0) {
	fprintf(fp, "Filesystem check: log file directory %s is writeable.\n", dname); fflush(fp);
      }
    }
    free(dirc);
  }

  if (gv->SEISMO > 0) {
    /* check seismogram directory */
    char *dirc = strdup(gv->SEIS_FILE);
    if (!dirc) declare_error("Could not copy string in check_fs.c - memory issue encountered.");
    char *dname = dirname(dirc);
    if (access(dname, W_OK) != 0) {
      fprintf(fp, "\n==================================================================\n");
      fprintf(fp, "ERROR: PE=%d, cannot write to seismogram directory %s.\n", MYID, dname);
      fprintf(fp, "\n==================================================================\n");
      fserr = 1;
    } else {
      if (MYID == 0) {
	fprintf(fp, "Filesystem check: seismogram directory %s is writeable.\n", dname); fflush(fp);
      }
    }
    free(dirc);
  }

  if (gv->SNAP > 0) {
    /* check snapshot directory */
    char *dirc = strdup(gv->SNAP_FILE);
    if (!dirc) declare_error("Could not copy string in check_fs.c - memory issue encountered.");
    char *dname = dirname(dirc);
    if (access(dname, W_OK) != 0) {
      fprintf(fp, "\n==================================================================\n");
      fprintf(fp, "ERROR: PE=%d, cannot write to snapshot directory %s.\n", MYID, dname);
      fprintf(fp, "\n==================================================================\n");
      fserr = 1;
    } else {
      if (MYID == 0) {
	fprintf(fp, "Filesystem check: snapshot directory %s is writeable.\n", dname); fflush(fp);
      }
    }
    free(dirc);
  }

  if (gv->WRITE_MODELFILES > 0) {
    /* check model file directory */
    char *dirc = strdup(gv->MFILE);
    if (!dirc) declare_error("Could not copy string in check_fs.c - memory issue encountered.");
    char *dname = dirname(dirc);
    if (access(dname, W_OK) != 0) {
      fprintf(fp, "\n==================================================================\n");
      fprintf(fp, "ERROR: PE=%d, cannot write to model file directory %s.\n", MYID, dname);
      fprintf(fp, "\n==================================================================\n");
      fserr = 1;
    } else {
      if (MYID == 0) {
	fprintf(fp, "Filesystem check: model file directory %s is writeable.\n", dname); fflush(fp);
      }
    }
    free(dirc);
  }

  if (1 == gv->SIGOUT) {
    /* check signal output directory */
    char *dirc = strdup(gv->SIGOUT_FILE);
    if (!dirc) declare_error("Could not copy string in check_fs.c - memory issue encountered.");
    char *dname = dirname(dirc);
    if (access(dname, W_OK) != 0) {
      fprintf(fp, "\n==================================================================\n");
      fprintf(fp, "ERROR: PE=%d, cannot write to source output directory %s.\n", MYID, dname);
      fprintf(fp, "\n==================================================================\n");
      fserr = 1;
    } else {
      if (MYID == 0) {
	fprintf(fp, "Filesystem check: source output directory %s is writeable.\n", dname); fflush(fp);
      }
    }
    free(dirc);
  }
  
  
  /********************************************/
  /* ERROR                                    */
  /********************************************/
  if (fserr) {
    fprintf(fp, "\n");
    sprintf(errmsg, "\n  in: <check_fs.c> \n");
    declare_error(errmsg);
  }
}
