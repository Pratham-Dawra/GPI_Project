
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

/* -------------------------------------------------------------
 * This is function freemem.
 * De-allocating memory of variables.
 *
 * -------------------------------------------------------------*/

#include "fd.h"

void freemem(MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{

    /* de-allocate buffer for messages */
    MPI_Buffer_detach(gv->BUFF_ADDR, &gv->BUFFSIZE);

    /* de-allocate buffer for subgrid arrays */
    if (gv->GY)
        free_ivector(gv->GY, 1, 4);
    if (gv->GX)
        free_ivector(gv->GX, 1, 4);

    /* de-allocate wavefield buffers */
    freemem_wavefield(mpw, gv);

    /* de-allocate model buffers */
    freemem_model(mpm, gv);

    /* de-allocate buffer for seismogram output, merged seismogram section of all PEs */
    if (gv->SEISMO_FULLDATA)
        free_matrix(gv->SEISMO_FULLDATA, 1, gv->NTRG, 1, gv->NS);

    if (gv->SECTIONVX)
        free_matrix(gv->SECTIONVX, 1, gv->NTR, 1, gv->NS);
    if (gv->SECTIONVY)
        free_matrix(gv->SECTIONVY, 1, gv->NTR, 1, gv->NS);
    if (gv->SECTIONP)
        free_matrix(gv->SECTIONP, 1, gv->NTR, 1, gv->NS);
    if (gv->SECTIONDIV)
        free_matrix(gv->SECTIONDIV, 1, gv->NTR, 1, gv->NS);
    if (gv->SECTIONCURL)
        free_matrix(gv->SECTIONCURL, 1, gv->NTR, 1, gv->NS);

}
