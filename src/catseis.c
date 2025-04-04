
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

/*-------------------------------------------------------------
 * Cat seismograms for collecitve output of seismograms
 *------------------------------------------------------------- */

#include "fd.h"
#include "util.h"

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns)
{
    /* temporary global data array for MPI-exchange */
    float **fulldata2 = matrix(1, ntr_glob, 1, ns);

    int k = 0;                  /* trace counter for local data array */

    /* loop over global traces: copy traces of local array  */
    /* to appropriate locations in the global array     */
    for (int i = 1; i <= ntr_glob; i++) {
        if (recswitch[i]) {
            k++;
            for (int j = 1; j <= ns; j++)
                fulldata2[i][j] = data[k][j];
        }
    }

    MPI_Allreduce(&fulldata2[1][1], &fulldata[1][1], ntr_glob * ns, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    free_matrix(fulldata2, 1, ntr_glob, 1, ns);
}
