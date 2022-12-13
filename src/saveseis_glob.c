
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
 *   write merged seismograms from all PEs to files,
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void saveseis_glob(float **sectiondata, int **recpos, int ntr,
                   float **srcpos, int ishot, int ns, int sectiondatatype, GlobVar *gv)
{
    char vxf[STRING_SIZE * 2], vyf[STRING_SIZE * 2], curlf[STRING_SIZE * 2], divf[STRING_SIZE * 2], pf[STRING_SIZE * 2];
    char file_ext[5];

    //switch (gv->SEIS_FORMAT[0]){
    switch (gv->SEIS_FORMAT) {
      case 0:
          sprintf(file_ext, "sgy");
          break;
      case 1:
          sprintf(file_ext, "su");
          break;
      case 2:
          sprintf(file_ext, "txt");
          break;
      case 3:
          sprintf(file_ext, "bin");
          break;
      case 4:
          sprintf(file_ext, "sgy");
          break;
      case 5:
          sprintf(file_ext, "sgy");
          break;
    }

    if (gv->RUN_MULTIPLE_SHOTS) {
        sprintf(vxf, "%s_vx.%s.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(vyf, "%s_vy.%s.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(curlf, "%s_curl.%s.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(divf, "%s_div.%s.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(pf, "%s_p.%s.shot%d", gv->SEIS_FILE, file_ext, ishot);
    } else {
        sprintf(vxf, "%s_vx.%s", gv->SEIS_FILE, file_ext);
        sprintf(vyf, "%s_vy.%s", gv->SEIS_FILE, file_ext);
        sprintf(curlf, "%s_curl.%s", gv->SEIS_FILE, file_ext);
        sprintf(divf, "%s_div.%s", gv->SEIS_FILE, file_ext);
        sprintf(pf, "%s_p.%s", gv->SEIS_FILE, file_ext);
    }

    switch (sectiondatatype) {
      case 1:                  /* particle velocities vx only */
          log_info("Writing %d merged seismogram traces (vx) to %s.\n", ntr, vxf);
          outseis_glob(fopen(vxf, "w"), sectiondata, recpos, ntr, srcpos, ns, gv->SEIS_FORMAT, ishot, 1, gv);
          break;
      case 2:                  /* particle velocities vy only */
          log_info("Writing %d merged seismogram traces (vy) to %s.\n", ntr, vyf);
          outseis_glob(fopen(vyf, "w"), sectiondata, recpos, ntr, srcpos, ns, gv->SEIS_FORMAT, ishot, 2, gv);
          break;
      case 4:                  /* pressure only */
          log_info("Writing %d merged seismogram traces (p) to %s.\n", ntr, pf);
          outseis_glob(fopen(pf, "w"), sectiondata, recpos, ntr, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 5:                  /* curl only */
          log_info("Writing %d merged seismogram traces (div) to %s.\n", ntr, divf);
          outseis_glob(fopen(divf, "w"), sectiondata, recpos, ntr, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 6:                  /* div only */
          log_info("Writing %d merged seismogram traces (curl) to %s.\n", ntr, curlf);
          outseis_glob(fopen(curlf, "w"), sectiondata, recpos, ntr, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
    }

    return;
}
