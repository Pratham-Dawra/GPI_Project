/*------------------------------------------------------------------------
 * Copyright (C) 2024 For the list of authors, see file AUTHORS.
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
 * along with IFOS. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 *------------------------------------------------------------------------*/

#ifndef SEISMO_SHIFT_H
#define SEISMO_SHIFT_H

#include "globvar_struct.h"
#include "fd.h"

void seismo_shift(st_seismogram *section,  // struct containing seismogram buffers
		  int ntr,                 // number of traces
		  float fc,                // main frequency of source (only used if gv->AUTO_SHIFT==1)
		  const GlobVar *gv);      // struct containing global parameters

#endif
