
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
 *   Write merged seismograms from all PEs to disk
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "su_struct.h"
#include "write_su.h"
#include "logging.h"

void outseis_glob(FILE * fpdata, float **section,
                  int **recpos, int ntr, float **srcpos, int ns, int seis_form, int ishot, int comp, GlobVar * gv)
{
    const float xshift = 800.0, yshift = 800.0;

    SUhead tr;
    float xr, yr, x, y;
    float XS = 0.0, YS = 0.0;

    switch (seis_form) {
      case 1:
          /* if there are more than one source position/coordinate specified
           * in SOURCE_FILE, only the first position is written into trace header fields,
           * in case of RUN_MULTIPLE_SHOTS is activated the individual shot position is used 
           */
          XS = srcpos[1][ishot];
          YS = srcpos[2][ishot];
          for (int tracl1 = 1; tracl1 <= ntr; tracl1++) {   /*SEGY (without file-header) */
              init_SUhead(&tr); /* set all headers to zero by default */
              xr = recpos[1][tracl1] * gv->DH;
              yr = recpos[2][tracl1] * gv->DH;
              x = xr - gv->REFREC[1];
              y = yr - gv->REFREC[2];

              tr.tracl = (int)recpos[3][tracl1];

              tr.ep = comp;
              tr.cdp = (int)recpos[3][tracl1];
              tr.trid = (short)1;   /* trace identification code: 1=seismic */
              tr.offset = (int)iround(sqrt((XS - xr) * (XS - xr) + (YS - yr) * (YS - yr)));
              tr.gelev = (int)iround(yr * 1000.0);
              tr.sdepth = (int)iround(YS * 1000.0); /* source depth (positive) */
              /* angle between receiver position and reference point
               * (sperical coordinate system: swdep=theta, gwdep=phi) */
              tr.swdep = iround(((360.0 / (2.0 * PI)) * atan2(x - xshift, y - yshift)) * 1000.0);
              tr.scalel = (short)-1000;
              tr.scalco = (short)-1000;
              tr.sx = (int)iround(XS * 1000.0); /* X source coordinate */

              /* group coordinates */
              tr.gx = (int)iround(xr * 1000.0);

              tr.ns = (unsigned short)ns;   /* number of samples in this trace */
              tr.dt = (unsigned short)iround(((float)gv->NDT * gv->DT) * 1.0e6);    /* sample interval in micro-seconds */
              tr.d1 = (float)tr.dt * 1.0e-6;    /* sample spacing for non-seismic data */

              tr.tracr = tr.tracl;  /* trace sequence number within reel */
              tr.fldr = ishot;  /* field record number */
              tr.tracf = tracl1;    /* trace number within field record */
              tr.ep = ishot;    /* energy source point number */

              // section is 1-based rather than zero-based
              su_write_trace(fpdata, &tr, &(section[tracl1][1]));
          }
          break;

      case 2:
          for (int i = 1; i <= ntr; i++) {  /* ASCII ONE COLUMN */
              for (int j = 1; j <= ns; j++)
                  fprintf(fpdata, "%e\n", section[i][j]);
          }
          break;

      case 3:                  /* BINARY */
          for (int i = 1; i <= ntr; i++)
              for (int j = 1; j <= ns; j++) {
                  fwrite(&section[i][j], sizeof(float), 1, fpdata);
              }
          break;

      default:
          log_warn("Unknown data format for seismograms. No output written!\n");
    }

    fclose(fpdata);
}
