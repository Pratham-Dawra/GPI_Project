
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the absorbing frame (ABS=2 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   GX and GY are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void update_s_acoustic_abs(int nt, MemModel * mpm, MemWavefield * mpw, GlobVar * gv)
{
    double time1 = 0.0, time2 = 0.0;

    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time1 = MPI_Wtime();
        log_debug("Updating stress components absorbing frame...\n");
    }
 
    if (gv->WEQ >= EL_ISO && gv->WEQ <= VEL_TTI) {  /* elastic cases */
        
        /* left boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* right boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* top boundary */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* bottom boundary */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* corners */
        
        /*left-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /*left-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* right-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
        
        /* right-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s(i, j, mpm, mpw);
            }
        }
    }
    else if (gv->WEQ == AC_VTI || gv->WEQ == AC_TTI)
    {
        /* left boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* right boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* top boundary */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* bottom boundary */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* corners */
        
        /*left-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /*left-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* right-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
        
        /* right-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac_vti(i, j, mpm, mpw);
                abs_update_s_ac1(i, j, mpm, mpw);
            }
        }
    }
    else if (gv->WEQ == AC_ISO)
    {
        /* left boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* right boundary */
        for (int j = gv->GY[2] + 1; j <= gv->GY[3]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* top boundary */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* bottom boundary */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[2] + 1; i <= gv->GX[3]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* corners */
        
        /*left-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /*left-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[1]; i <= gv->GX[2]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* right-top */
        for (int j = gv->GY[1]; j <= gv->GY[2]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
        
        /* right-bottom */
        for (int j = gv->GY[3] + 1; j <= gv->GY[4]; j++) {
            for (int i = gv->GX[3] + 1; i <= gv->GX[4]; i++) {
                wavefield_update_s_ac(i, j, mpm, mpw);
                abs_update_s_ac2(i, j, mpm, mpw);
            }
        }
    }
    if ((gv->MPID == 0) && ((nt + (gv->OUTNTIMESTEPINFO - 1)) % gv->OUTNTIMESTEPINFO) == 0) {
        time2 = MPI_Wtime();
        log_debug("Finished updating stress components (real time: %4.3fs).\n", time2 - time1);
    }
}
