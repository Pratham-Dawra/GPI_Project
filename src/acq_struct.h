
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
 * globvar_struct.h - global variables of viscoelastic 2D FD program
 * generally, for the names of the global variables uppercase letters are used
 *
 * ----------------------------------------------------------------------*/

#ifndef ACQ_STRUCT_H_INCLUDED
#define ACQ_STRUCT_H_INCLUDED

typedef struct {
    // Acquisition parameters
    float **srcpos;
    float **srcpos_loc;
    int *srcswitch;

    int **recpos;
    int **recpos_loc;
    int *recswitch;

    int nsrc;                   // number of shots (in case of multiple shots)
    int nsrc_loc;               // number of shots (for each subdomain)
    float **signals;            // source wavelet

} AcqVar;

#endif
