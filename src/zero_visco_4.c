
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

void zero_visco_4(int nx1, int nx2, int ny1, int ny2, float **vxx_1, float **vxx_2, float **vxx_3, float **vxx_4,
                  float **vyy_1, float **vyy_2, float **vyy_3, float **vyy_4, float **vxy_1, float **vxy_2,
                  float **vxy_3, float **vxy_4, float **vyx_1, float **vyx_2, float **vyx_3, float **vyx_4,
                  float **svx_1, float **svx_2, float **svx_3, float **svx_4, float **svy_1, float **svy_2,
                  float **svy_3, float **svy_4, float ***pr_2, float ***pr_3, float ***pr_4, float ***pp_2,
                  float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4, GlobVar *gv)
{
    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            vxx_1[j][i] = 0.0f;
            vxx_2[j][i] = 0.0f;
            vxx_3[j][i] = 0.0f;
            vxx_4[j][i] = 0.0f;
            vyy_1[j][i] = 0.0f;
            vyy_2[j][i] = 0.0f;
            vyy_3[j][i] = 0.0f;
            vyy_4[j][i] = 0.0f;
            vxy_1[j][i] = 0.0f;
            vxy_2[j][i] = 0.0f;
            vxy_3[j][i] = 0.0f;
            vxy_4[j][i] = 0.0f;
            vyx_1[j][i] = 0.0f;
            vyx_2[j][i] = 0.0f;
            vyx_3[j][i] = 0.0f;
            vyx_4[j][i] = 0.0f;
            svx_1[j][i] = 0.0f;
            svx_2[j][i] = 0.0f;
            svx_3[j][i] = 0.0f;
            svx_4[j][i] = 0.0f;
            svy_1[j][i] = 0.0f;
            svy_2[j][i] = 0.0f;
            svy_3[j][i] = 0.0f;
            svy_4[j][i] = 0.0f;
        }
    }

    for (int j = ny1; j <= ny2; j++) {
        for (int i = nx1; i <= nx2; i++) {
            for (int l = 1; l <= gv->L; l++) {
                pr_2[j][i][l] = 0.0f;
                pr_3[j][i][l] = 0.0f;
                pr_4[j][i][l] = 0.0f;

                pp_2[j][i][l] = 0.0f;
                pp_3[j][i][l] = 0.0f;
                pp_4[j][i][l] = 0.0f;

                pq_2[j][i][l] = 0.0f;
                pq_3[j][i][l] = 0.0f;
                pq_4[j][i][l] = 0.0f;
            }
        }
    }
}
