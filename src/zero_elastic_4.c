
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
#include "fd.h"

void zero_elastic_4(int nx1, int nx2, int ny1, int ny2, MemWavefield * mpw)
{
    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            mpw->vxx_1[j][i] = 0.0f;
            mpw->vxx_2[j][i] = 0.0f;
            mpw->vxx_3[j][i] = 0.0f;
            mpw->vxx_4[j][i] = 0.0f;
            mpw->vyy_1[j][i] = 0.0f;
            mpw->vyy_2[j][i] = 0.0f;
            mpw->vyy_3[j][i] = 0.0f;
            mpw->vyy_4[j][i] = 0.0f;
            mpw->vxy_1[j][i] = 0.0f;
            mpw->vxy_2[j][i] = 0.0f;
            mpw->vxy_3[j][i] = 0.0f;
            mpw->vxy_4[j][i] = 0.0f;
            mpw->vyx_1[j][i] = 0.0f;
            mpw->vyx_2[j][i] = 0.0f;
            mpw->vyx_3[j][i] = 0.0f;
            mpw->vyx_4[j][i] = 0.0f;
            mpw->svx_1[j][i] = 0.0f;
            mpw->svx_2[j][i] = 0.0f;
            mpw->svx_3[j][i] = 0.0f;
            mpw->svx_4[j][i] = 0.0f;
            mpw->svy_1[j][i] = 0.0f;
            mpw->svy_2[j][i] = 0.0f;
            mpw->svy_3[j][i] = 0.0f;
            mpw->svy_4[j][i] = 0.0f;
        }
    }
}
