
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
 * perform_struct.h - variables to measure performance of viscoelastic 2D
 * FD program
 * ----------------------------------------------------------------------*/

#ifndef PERFORM_STRUCT_H_INCLUDED
#define PERFORM_STRUCT_H_INCLUDED

typedef struct {
    int infocounter;
    double time_av_v_update;
    double time_av_s_update;
    double time_av_v_exchange;
    double time_av_s_exchange;
    double time_av_timestep;
} Perform;

#endif
