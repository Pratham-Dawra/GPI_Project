
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

void saveseis_glob(float **sectiondata, int **recpos, float **srcpos, int ishot, int ns, int sectiondatatype,
                   GlobVar *gv)
{
    char vxf[STRING_SIZE * 2], vyf[STRING_SIZE * 2], curlf[STRING_SIZE * 2], divf[STRING_SIZE * 2], pf[STRING_SIZE * 2];
    char vxdifff[STRING_SIZE * 2], vydifff[STRING_SIZE * 2], pdifff[STRING_SIZE * 2];
    char vxobsf[STRING_SIZE * 2], vyobsf[STRING_SIZE * 2], pobsf[STRING_SIZE * 2];
    char vxsynf[STRING_SIZE * 2], vysynf[STRING_SIZE * 2], psynf[STRING_SIZE * 2];
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
    if (gv->MODE == FWI) {
        sprintf(vxdifff, "%s_vx.%s.diff.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(vydifff, "%s_vy.%s.diff.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(pdifff, "%s_p.%s.diff.shot%d", gv->SEIS_FILE, file_ext, ishot);

        sprintf(vxobsf, "%s_vx.%s.obs.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(vyobsf, "%s_vy.%s.obs.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(pobsf, "%s_p.%s.obs.shot%d", gv->SEIS_FILE, file_ext, ishot);

        sprintf(vxsynf, "%s_vx.%s.syn.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(vysynf, "%s_vy.%s.syn.shot%d", gv->SEIS_FILE, file_ext, ishot);
        sprintf(psynf, "%s_p.%s.syn.shot%d", gv->SEIS_FILE, file_ext, ishot);
    }

    switch (sectiondatatype) {
      case 1:                  /* particle velocities vx only */
          log_info("Writing %d merged seismogram traces (vx) to %s.\n", gv->NTRG, vxf);
          outseis_glob(fopen(vxf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 1, gv);
          break;
      case 2:                  /* particle velocities vy only */
          log_info("Writing %d merged seismogram traces (vy) to %s.\n", gv->NTRG, vyf);
          outseis_glob(fopen(vyf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 2, gv);
          break;
      case 4:                  /* pressure only */
          log_info("Writing %d merged seismogram traces (p) to %s.\n", gv->NTRG, pf);
          outseis_glob(fopen(pf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 5:                  /* curl only */
          log_info("Writing %d merged seismogram traces (div) to %s.\n", gv->NTRG, divf);
          outseis_glob(fopen(divf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 6:                  /* div only */
          log_info("Writing %d merged seismogram traces (curl) to %s.\n", gv->NTRG, curlf);
          outseis_glob(fopen(curlf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 11:                 /* difference of particle velocities vx only */
          log_info("Writing %d merged seismogram traces (vx) to %s.\n", gv->NTRG, vxdifff);
          outseis_glob(fopen(vxdifff, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 1, gv);
          break;
      case 12:                 /* difference of particle velocities vy only */
          log_info("Writing %d merged seismogram traces (vy) to %s.\n", gv->NTRG, vydifff);
          outseis_glob(fopen(vydifff, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 2, gv);
          break;
      case 14:                 /* difference of pressure only */
          log_info("Writing %d merged seismogram traces (p) to %s.\n", gv->NTRG, pdifff);
          outseis_glob(fopen(pdifff, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 21:                 /* particle velocities vx of filtered field data */
          log_info("Writing %d merged seismogram traces (vx) to %s.\n", gv->NTRG, vxobsf);
          outseis_glob(fopen(vxobsf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 1, gv);
          break;
      case 22:                 /* particle velocities vy of filtered field data */
          log_info("Writing %d merged seismogram traces (vy) to %s.\n", gv->NTRG, vyobsf);
          outseis_glob(fopen(vyobsf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 2, gv);
          break;
      case 24:                 /* pressure of filtered field data */
          log_info("Writing %d merged seismogram traces (p) to %s.\n", gv->NTRG, pobsf);
          outseis_glob(fopen(pobsf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
      case 31:                 /* particle velocities vx of filtered synthetic data */
          log_info("Writing %d merged seismogram traces (vx) to %s.\n", gv->NTRG, vxsynf);
          outseis_glob(fopen(vxsynf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 1, gv);
          break;
      case 32:                 /* particle velocities vy of filtered synthetic data */
          log_info("Writing %d merged seismogram traces (vy) to %s.\n", gv->NTRG, vysynf);
          outseis_glob(fopen(vysynf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 2, gv);
          break;
      case 34:                 /* pressure of filtered synthetic data */
          log_info("Writing %d merged seismogram traces (p) to %s.\n", gv->NTRG, psynf);
          outseis_glob(fopen(psynf, "w"), sectiondata, recpos, gv->NTRG, srcpos, ns, gv->SEIS_FORMAT, ishot, 0, gv);
          break;
    }

    return;
}
