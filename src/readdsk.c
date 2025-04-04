
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
 *   Read one single amplitude from file                                   
 *  ----------------------------------------------------------------------*/

#include <stdio.h>
#include "fd.h"
#include "logging.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/

float readdsk(FILE * fp_in, int format)
{
    float amp = 0.0f;

    switch (format) {
      case 1:                  /* SU */
          log_fatal("Sorry, SU-format for snapshots not implemented yet.\n");
          break;
      case 2:
          /*ASCII*/ fscanf(fp_in, "%e\n", &amp);
          break;
      case 3:                  /* BINARY */
          fread(&amp, sizeof(float), 1, fp_in);
          break;
      default:
          log_error("Unknown format for the snapshot data.\n");
          log_fatal("Unsupported format, no files processed.\n");
    }
    return amp;
}
