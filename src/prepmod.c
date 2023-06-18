
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
 *   Prepare model parameters
 * ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void prepmod(MemModel *mpm, GlobVar *gv)
{
    /* ------------------------------------------------------------------------
     *   For the calculation of the material parameters between gridpoints
     *   they have to be averaged. For this, values lying at 0 and NX+1,
     *   for example, are required on the local grid. These are now copied from the
     *   neighbouring grids
     * ----------------------------------------------------------------------*/

    switch (gv->WEQ) {
      case AC_ISO:             /* acoustic */
          matcopy_elastic(mpm->prho, mpm->ppi, mpm->ppi, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
          break;
      case EL_ISO:             /* elastic */
          matcopy_elastic(mpm->prho, mpm->ppi, mpm->pu, gv);
          av_mue(mpm->pu, mpm->puipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
          break;
      case VEL_ISO:            /* viscoelastic */
          matcopy(mpm->prho, mpm->ppi, mpm->pu, mpm->ptaus, mpm->ptaup, gv);
          av_mue(mpm->pu, mpm->puipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
          av_tau(mpm->ptaus, mpm->ptausipjp, gv);
          break;
      case EL_VTI:             /* elastic VTI */
          matcopy_elastic(mpm->prho, mpm->pc11, mpm->pc55, gv);
          av_mue(mpm->pc55, mpm->pc55ipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
#ifdef EBUG
          debug_check_matrix(mpm->pc55ipjp, 0, gv->NX, gv->NY, 55, 0, "pc55ipjp");
          debug_check_matrix(mpm->prip, 0, gv->NX, gv->NY, 55, 0, "prip");
          debug_check_matrix(mpm->prjp, 0, gv->NX, gv->NY, 55, 0, "prjp");
#endif
          break;
      case VEL_VTI:            /* viscoelastic VTI */
          matcopy_elastic(mpm->prho, mpm->ptau55, mpm->pc55, gv);
          av_mue(mpm->pc55, mpm->pc55ipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
          av_tau(mpm->ptau55, mpm->ptau55ipjp, gv);
#ifdef EBUG
          debug_check_matrix(mpm->pc55ipjp, 0, gv->NX, gv->NY, 66, 0, "pc55ipjp");
          debug_check_matrix(mpm->prip, 0, gv->NX, gv->NY, 66, 0, "prip");
          debug_check_matrix(mpm->prjp, 0, gv->NX, gv->NY, 66, 0, "prjp");
          debug_check_matrix(mpm->ptau55ipjp, 0, gv->NX, gv->NY, 66, 0, "ptau55ipjp");
#endif
          break;
      case EL_TTI:             /* elastic TTI */
          matcopy_elastic(mpm->prho, mpm->pc11, mpm->pc55, gv);
          matcopy_elastic(mpm->prho, mpm->pc15, mpm->pc35, gv);
          av_mue(mpm->pc55, mpm->pc55ipjp, gv);
          av_mue(mpm->pc15, mpm->pc15ipjp, gv);
          av_mue(mpm->pc35, mpm->pc35ipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
#ifdef EBUG
          debug_check_matrix(mpm->pc55ipjp, 0, gv->NX, gv->NY, 77, 0, "pc55ipjp");
          debug_check_matrix(mpm->prip, 0, gv->NX, gv->NY, 77, 0, "prip");
          debug_check_matrix(mpm->prjp, 0, gv->NX, gv->NY, 77, 0, "prjp");
          debug_check_matrix(mpm->pc15ipjp, 0, gv->NX, gv->NY, 77, 0, "pc15ipjp");
          debug_check_matrix(mpm->pc35ipjp, 0, gv->NX, gv->NY, 77, 0, "pc35ipjp");
#endif
          break;
      case VEL_TTI:            /* viscoelastic TTI */
          matcopy_elastic(mpm->prho, mpm->ptau55, mpm->pc55, gv);
          matcopy_elastic(mpm->prho, mpm->ptau15, mpm->pc15, gv);
          matcopy_elastic(mpm->prho, mpm->ptau35, mpm->pc35, gv);
          av_mue(mpm->pc55, mpm->pc55ipjp, gv);
          av_mue(mpm->pc15, mpm->pc15ipjp, gv);
          av_mue(mpm->pc35, mpm->pc35ipjp, gv);
          av_rho(mpm->prho, mpm->prip, mpm->prjp, gv);
          av_tau(mpm->ptau55, mpm->ptau55ipjp, gv);
          av_tau(mpm->ptau15, mpm->ptau15ipjp, gv);
          av_tau(mpm->ptau35, mpm->ptau35ipjp, gv);
#ifdef EBUG
          debug_check_matrix(mpm->pc55ipjp, 0, gv->NX, gv->NY, 88, 0, "pc55ipjp");
          debug_check_matrix(mpm->prip, 0, gv->NX, gv->NY, 88, 0, "prip");
          debug_check_matrix(mpm->prjp, 0, gv->NX, gv->NY, 88, 0, "prjp");
          debug_check_matrix(mpm->pc15ipjp, 0, gv->NX, gv->NY, 88, 0, "pc15ipjp");
          debug_check_matrix(mpm->pc35ipjp, 0, gv->NX, gv->NY, 88, 0, "pc35ipjp");
          debug_check_matrix(mpm->ptau55ipjp, 0, gv->NX, gv->NY, 88, 0, "ptau55ipjp");
          debug_check_matrix(mpm->ptau15ipjp, 0, gv->NX, gv->NY, 88, 0, "ptau15ipjp");
          debug_check_matrix(mpm->ptau35ipjp, 0, gv->NX, gv->NY, 88, 0, "ptau35ipjp");
#endif
          break;
      case VAC_ISO:            /* viscoacoustic */
          log_fatal("not yet implemented\n");
          break;
      default:
          log_fatal("Unknown WEQ.\n");
    }

    /* Preparing memory variables for update_s (viscoelastic only) */

    if (gv->FDORDER_TIME == 2) {
        switch (gv->WEQ) {
          case AC_ISO:         /* acoustic */
              prepare_update_s_ac(mpm, gv);
              break;
          case EL_ISO:         /* elastic */
              prepare_update_s_el(mpm, gv);
              break;
          case VEL_ISO:        /* viscoelastic */
              prepare_update_s_visc(mpm, gv);
              break;
          case EL_VTI:
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc11);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc33);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc13);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc55ipjp);
              break;
          case VEL_VTI:        /* viscoelastic VTI */
              prepare_update_s_vti(mpm, gv);
              break;
          case EL_TTI:
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc11);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc33);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc13);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc15);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc35);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc55ipjp);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc35ipjp);
              dt_mult(gv->NX, gv->NY, gv->DT, mpm->pc15ipjp);
              break;
          case VEL_TTI:        /* viscoelastic TTI */
              prepare_update_s_tti(mpm, gv);
              break;
          case VAC_ISO:        /* viscoacoustic */
              log_fatal("not yet implemented\n");
              break;
          default:
              log_fatal("Unknown WEQ.\n");
        }
    }

    gv->DTDH = gv->DT / gv->DH;

    if ((gv->WEQ == VEL_ISO) && gv->FDORDER_TIME == 4) {
        prepare_update_s_4(mpm, gv);
    }
}
