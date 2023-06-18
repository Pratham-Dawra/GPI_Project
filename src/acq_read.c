
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 * Allocate acquistion parameters
 *----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

int acq_read(AcqVar *acq, GlobVar *gv)
{
    int nshots = 1;

    int *topo = NULL;
    if (1 == gv->SOURCE_TOPO || 1 == gv->REC_TOPO) {
	/* Scan for topography in P-wave velocity model */
	if (0 == gv->MPID) {
	    log_info("------------------------- Topography scan -------------------\n");
	    topo = scan_topo(gv);
	    int min_topo = 999999;
	    int min_tidx = 1;
	    int max_topo = 0;
	    int max_tidx = 1;
	    for (int i=1; i<=gv->NXG; ++i) {
		if (topo[i] < min_topo) {
		    min_topo = topo[i];
		    min_tidx = i;
		}
		if (topo[i] > max_topo) {
		    max_topo = topo[i];
		    max_tidx = i;
		}
	    }
	    log_info("Maximum topography at x=%.4fm, %.4fm below top of model.\n",
		     (min_tidx-1)*gv->DH, (min_topo-1)*gv->DH);
	    log_info("Minimum topography at x=%.4fm, %.4fm below top of model.\n",
		     (max_tidx-1)*gv->DH, (max_topo-1)*gv->DH);
	} else {
	    topo = ivector(1, gv->NXG);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&(topo[1]), gv->NXG, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Reading receiver positions from REC_FILE */
    if (gv->SEISMO) {
        acq->recpos = receiver(gv, topo);
        acq->recswitch = ivector(1, gv->NTRG);
        acq->recpos_loc = splitrec(acq->recpos, acq->recswitch, gv);
    }

    /* Reading source positions from SOURCE_FILE */
    sources(acq, gv, topo);

    if (gv->RUN_MULTIPLE_SHOTS) {
        nshots = acq->nsrc;
    }

    return nshots;
}
