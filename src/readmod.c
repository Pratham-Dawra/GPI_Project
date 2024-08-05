
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

/* ------------------------------------------------------------------------
 * Read model properties from files
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void readmod(MemModel *mpm, MemInv *minv, GlobVar *gv, GlobVarInv *vinv)
{
    /* create model grids */
    if (gv->READMOD) {
        switch (gv->WEQ) {
          case AC_ISO:         /* acoustic */
              readmod_acoustic(mpm, gv);
              break;
          case EL_ISO:         /* elastic */
              readmod_elastic(mpm, minv, gv);
              break;
          case VEL_ISO:        /* viscoelastic */
              readmod_visco(mpm, minv, gv);
              break;
          case EL_VTI:         /* elastic VTI */
              readmod_elastic_vti(mpm, minv, gv);
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              readmod_visco_vti(mpm, minv, gv);
              break;
          case EL_TTI:         /* elastic TTI */
              readmod_elastic_tti(mpm, minv, gv);
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              readmod_visco_tti(mpm, minv, gv);
              break;
          case VAC_ISO:        /* viscoacoustic */
              log_fatal("not yet implemented\n");
              break;
          default:
              log_fatal("Unknown WEQ.\n");
        }
    } else {
        switch (gv->WEQ) {
          case EL_ISO:         /* elastic */
              model_elastic(mpm, gv);
              break;
          case VEL_ISO:        /* viscoelastic */
              model_visco(mpm, gv);
              break;
          case EL_VTI:         /* elastic VTI */
              model_elastic_VTI(mpm, gv);
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              model_visco_vti(mpm, gv);
              break;
          case EL_TTI:         /* elastic TTI */
              model_elastic_TTI(mpm, gv);
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              model_visco_tti(mpm, gv);
              break;
          default:
              log_fatal("Internal model for your chosen WEQ not implemented.\n");
              break;
        }
    }
    
    /* Some initial calculations for FWI */
    if (gv->MODE == FWI) {
        /* Get average values from material parameters */
        vinv->VP_AVG = average_matrix(minv->Vp0, gv);
        vinv->VS_AVG = average_matrix(minv->Vs0, gv);
        vinv->RHO_AVG = average_matrix(minv->Rho0, gv);

        //log_info("MYID = %d \t VP_AVG = %e \t VS_AVG = %e \t RHO_AVG = %e\n", gv->MPID, vinv->VP_AVG, vinv->VS_AVG, vinv->RHO_AVG);

        vinv->C_VP = vinv->VP_AVG * vinv->VP_AVG;
        vinv->C_VS = vinv->VS_AVG * vinv->VS_AVG;
        vinv->C_RHO = vinv->RHO_AVG * vinv->RHO_AVG;
    }
}
