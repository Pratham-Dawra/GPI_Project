
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
 *   Allocate acquistion parameters
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

int acq_read(AcqVar *acq, GlobVar *gv)
{
    int nshots = 1;

    /* Reading receiver positions from REC_FILE */
    if (gv->SEISMO) {
        acq->recpos = receiver(gv);
        acq->recswitch = ivector(1, gv->NTRG);
        acq->recpos_loc = splitrec(acq->recpos, acq->recswitch, gv);
    }

    /* Reading source positions from SOURCE_FILE */
    sources(acq, gv);

    if (gv->RUN_MULTIPLE_SHOTS) {
        nshots = acq->nsrc;
    }
    
    if (gv->SNAPSHOT_END == -9999 || gv->SNAPSHOT_END > nshots)
        gv->SNAPSHOT_END = nshots;

    return nshots;
}
