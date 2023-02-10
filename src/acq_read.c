
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
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "macros.h"

void acq_read(AcqVar *acq, GlobVar *gv){
  
    /* Reading receiver positions from REC_FILE */
    if (gv->SEISMO) {
        acq->recpos = receiver(gv);
        acq->recswitch = ivector(1, gv->NTRG);
        acq->recpos_loc = splitrec(acq->recpos, acq->recswitch, gv);
    }

    /* Memory allocation for saving the current source position */
    acq->srcpos_current = matrix(1, NSPAR, 1, 1);

    /* Reading source positions from SOURCE_FILE */
    sources(acq, gv);
  
}
