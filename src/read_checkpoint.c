
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

/* ------------------------------------------------------------------------
 * reads checkpoint file to continue a previously started modeling
 * (to allow for splitting of modeling jobs)
 * ------------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float **vx, float **vy, float **sxx, float **syy, float **sxy, GlobVar * gv)
{
    char myid[5];
    char checkptfile[STRING_SIZE];

    sprintf(checkptfile, "%s", gv->CHECKPTFILE);
    sprintf(myid, ".%d", gv->MPID);
    strcat(checkptfile, myid);

    FILE *fp = fopen(checkptfile, "rb");
    if (!fp)
        log_fatal("CHECKPTFILE %s cannot be opened for reading.\n", checkptfile);

    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            fread(&vx[j][i], sizeof(float), 1, fp);
            fread(&vy[j][i], sizeof(float), 1, fp);
            fread(&sxx[j][i], sizeof(float), 1, fp);
            fread(&syy[j][i], sizeof(float), 1, fp);
            fread(&sxy[j][i], sizeof(float), 1, fp);
        }
    }

    fclose(fp);

}
