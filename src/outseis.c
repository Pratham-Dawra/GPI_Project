
/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 *
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 *   Write seismograms to disk
 *----------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>
#include "fd.h"
#include "macros.h"
#include "su_struct.h"
#include "logging.h"
#include "write_su.h"

void outseis(FILE *fpdata, int comp, float **section,
             int **recpos, int **recpos_loc, int ntr, float **srcpos, int nsrc, GlobVar *gv)
{
    float xr, yr, zr, tfloat, xs = 0.0, ys = 0.0, zs = 0.0;
    SUhead tr;

    if (nsrc == 1) {
        /* only if one source position is specified in SOURCE_FILE,
         * source coordinates are written into trace header fields */
        xs = srcpos[1][1];
        ys = srcpos[2][1];
        zs = srcpos[3][1];
    }

    // number of actual output samples
    int ns = gv->NS;
    if (gv->NDT != 1) {
      ns = (int)ceilf((float)gv->NS / (float)gv->NDT);
    }

    const short scalefac = 1000;
    float *trace = NULL;

    switch (gv->SEIS_FORMAT) {
      case 0:
          /* FALLTHRU */
      case 1:                  // SU
	  if (gv->NDT != 1) {
	      trace = (float*)calloc(ns, sizeof(float));
	  }
          for (int tracl = 1; tracl <= ntr; tracl++) {
              init_SUhead(&tr); /* set all headers to zero by default */

	      if (recpos_loc) {
		  xr = (recpos[1][recpos_loc[4][tracl]]-1) * gv->DX;
		  yr = (recpos[2][recpos_loc[4][tracl]]-1) * gv->DY;
		  zr = (recpos[3][recpos_loc[4][tracl]]-1) * gv->DZ;
	      } else {
		  xr = (recpos[1][tracl]-1) * gv->DX;
		  yr = (recpos[2][tracl]-1) * gv->DY;
		  zr = (recpos[3][tracl]-1) * gv->DZ;
	      }
              float y = yr - gv->REFREC[2];
              float z = zr - gv->REFREC[3];

	      if (recpos_loc) {
		  tr.tracl = (int)recpos_loc[4][tracl];
	      } else {
		  tr.tracl = tracl;
	      }
              tr.ep = comp;
              tr.trid = (short)1;
              tr.offset = (int)iround(sqrt((xs-xr)*(xs-xr)+(ys-yr)*(ys-yr)+(zs-zr)*(zs-zr)));
              tr.gelev = iround(yr * scalefac);
              tr.sdepth = iround(ys * scalefac);
              // angle between receiver position and reference point
              // (spherical coordinate system: used for tunnel geometry)
              tr.gdel = iround(atan2(-y, z) * 180.0 * scalefac / PI);
              tr.gwdep = iround(sqrt(z * z + y * y) * scalefac);
              tr.scalel = (short)(-scalefac);
              tr.scalco = (short)(-scalefac);
              tr.sx = (int)iround(xs * scalefac);   // X source coordinate
              tr.sy = (int)iround(zs * scalefac);   // Z source coordinate
              tr.gx = iround(xr * scalefac);
              tr.gy = iround(zr * scalefac);
	      tr.ns = (unsigned short)ns;
              tr.dt = (unsigned short)iround(((float)gv->NDT * gv->DT) * 1.0e6);
              tr.d1 = (float)tr.dt * 1.0e-6;

	      if (gv->NDT != 1) {
		  int i = 0;
		  for (int j = 0; j <= gv->NS-1; j = j+gv->NDT) {
		      trace[i++] = section[tracl][j];
		  }
		  assert(i == ns); // cross-check
		  su_write_trace(fpdata, &tr, &(trace[0]));
	      } else {
		  su_write_trace(fpdata, &tr, &(section[tracl][0]));
	      }
          }
          break;
      case 2:                  // ASCII ONE COLUMN PER TRACE
          for (int j = 0; j <= gv->NS-1; j = j+gv->NDT) {
	      for (int i = 1; i <= ntr; i++) {
                  fprintf(fpdata, "%e\t", section[i][j]);
              }
              fprintf(fpdata, "\n");
          }
          break;
      case 3:                  // BINARY
          if (!gv->LITTLEBIG) {     // native endian
	      if (gv->NDT != 1) {
		  trace = (float*)calloc(ns, sizeof(float));
	      }
              for (int i = 1; i <= ntr; i++) {
	          if (gv->NDT != 1) {
		      int k = 0;
		      for (int j = 0; j <= gv->NS-1; j = j+gv->NDT) {
			  trace[k++] = section[i][j];
		      }
                      assert(k == ns); // cross-check
		      fwrite(&trace[0], sizeof(float), ns, fpdata);
		  } else {
		      fwrite(&section[i][0], sizeof(float), ns, fpdata);
		  }
	      }
          } else {              // foreign endian
              int *pint;
              for (int i = 1; i <= ntr; i++) {
                  for (int j = 0; j <= gv->NS-1; j = j+gv->NDT) {
                      tfloat = section[i][j];
                      pint = (int *)&tfloat;
                      *pint =
                          ((*pint >> 24) & 0xff) | ((*pint & 0xff) << 24) | ((*pint >> 8) & 0xff00) | ((*pint & 0xff00)
                                                                                                       << 8);
                      fwrite(&tfloat, sizeof(float), 1, fpdata);
                  }
              }
          }
          break;
      default:
          log_error("Unknown format for output seismograms. No output written.\n");
    }

    if (trace) {
        free(trace);
    }
    fclose(fpdata);
}
