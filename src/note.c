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

#include "fd.h"

void note(FILE *fp, GlobVar *gv) {

  int MYID;
  MPI_Comm_rank(MPI_COMM_WORLD, &MYID);

  int len = strlen(gv->LOG_FILE)-1;

  fprintf(fp,"\nPlease note:\n");
  fprintf(fp,"Each processing element (PE) is writing log information during program\n");
  fprintf(fp,"execution to %.*sPE.\n", len, gv->LOG_FILE);
  fprintf(fp,"See corresponding log files for further information on program status.\n");
  fprintf(fp,"Information about overall program execution\n");
  fprintf(fp,"(numerical artefacts, accuracy, computing times etc)\n");
  if (gv->LOG) {
    fprintf(fp,"will be written by PE 0 to");
    if ((gv->LOG == 1) || (gv->LOG == 3)) {
      fprintf(fp," standard output. \n");
    } else if (gv->LOG == 2) {
      fprintf(fp," %.*s%i.\n", len, gv->LOG_FILE, MYID);
    }
  } else {
    fprintf(fp," will NOT be output.");
  }
}
