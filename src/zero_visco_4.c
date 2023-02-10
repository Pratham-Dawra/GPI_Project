
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

void zero_visco_4(int j, int i, int l,MemWavefield * mpw)
{
    mpw->pr_2[j][i][l] = 0.0f;
    mpw->pr_3[j][i][l] = 0.0f;
    mpw->pr_4[j][i][l] = 0.0f;

    mpw->pp_2[j][i][l] = 0.0f;
    mpw->pp_3[j][i][l] = 0.0f;
    mpw->pp_4[j][i][l] = 0.0f;

    mpw->pq_2[j][i][l] = 0.0f;
    mpw->pq_3[j][i][l] = 0.0f;
    mpw->pq_4[j][i][l] = 0.0f;
}
