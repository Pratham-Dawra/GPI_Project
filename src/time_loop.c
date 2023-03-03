
/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
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
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   time loop or finite-difference forward modelling
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"

void time_loop(int ishot, int lsnap, int nsnap, float *hc, AcqVar *acq, MemModel *mpm,
               MemWavefield *mpw, GlobVar *gv, Perform *perf)
{
    int nt;
    double time3 = 0.0, time4 = 0.0, time5 = 0.0, time6 = 0.0, time7 = 0.0, time8 = 0.0;
    int lsamp = gv->NDT;

    for (nt = 1; nt <= gv->NT; nt++) {

        if (isnan(mpw->pvy[gv->NY / 2][gv->NX / 2])) {
            log_error("Time step: %d; pvy: %f.\n", nt, mpw->pvy[gv->NY / 2][gv->NX / 2]);
            log_fatal("Simulation is unstable!\n");
        }

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            log_info("Computing time step %d of %d.\n", nt, gv->NT);
            time3 = MPI_Wtime();
        }

        /*---------------------------------------------------------------
         * update of particle velocities --------------------------------
         *---------------------------------------------------------------*/
        if (gv->FDORDER_TIME == 2) {
            update_v_interior(nt, acq->srcpos_loc, acq->signals, acq->nsrc_loc, hc, mpm, mpw, gv);
#ifdef EBUG
            debug_check_matrix(mpw->pvx, nt, gv->NX, gv->NY, 121, 0, "pvx");
            debug_check_matrix(mpw->pvy, nt, gv->NX, gv->NY, 121, 0, "pvy");
#endif

            if (gv->FW) {
                if (gv->ABS_TYPE == 1) {
                    update_v_PML(gv->NX, gv->NY, nt, mpm, mpw, gv);
                }

                if (gv->ABS_TYPE == 2) {
                    update_v_abs(mpm, mpw, gv);
                }
#ifdef EBUG
                debug_check_matrix(mpw->pvx, nt, gv->NX, gv->NY, 122, 0, "pvx");
                debug_check_matrix(mpw->pvy, nt, gv->NX, gv->NY, 122, 0, "pvy");
#endif
            }
        }

        if (gv->FDORDER_TIME == 4) {
            update_v_interior_4(nt, acq->srcpos_loc, acq->signals, acq->nsrc_loc, hc, mpm, mpw, gv);
            if (gv->FW) {
                if (gv->ABS_TYPE == 1) {
                    update_v_PML_4(gv->NX, gv->NY, nt, mpm, mpw, gv);
                }

                if (gv->ABS_TYPE == 2) {
                    update_v_abs_4(nt, mpm, mpw, gv);
                }
            }

            /* Shift spatial derivations of the stress one time step back */
            shift_var2(&(mpw->svx_1), &(mpw->svx_2), &(mpw->svx_3), &(mpw->svx_4));
            shift_var2(&(mpw->svy_1), &(mpw->svy_2), &(mpw->svy_3), &(mpw->svy_4));
        }

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            time4 = MPI_Wtime();
            perf->time_av_v_update += (time4 - time3);
            log_debug("Starting particle velocity exchange between PEs...\n");
        }

        /*---------------------------------------------------------------*
         * ------- exchange of particle velocities between PEs --------- *
         *---------------------------------------------------------------*/
        exchange_v(mpw->pvx, mpw->pvy, mpw, gv);

        if ((gv->WEQ == EL_TTI) || (gv->WEQ == VEL_TTI)) {  /* TTI */
            v_derivatives(mpw, gv);
            exchange_v(mpw->pvxx, mpw->pvyy, mpw, gv);
            exchange_v(mpw->pvyx, mpw->pvxy, mpw, gv);
        }

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            time5 = MPI_Wtime();
            perf->time_av_v_exchange += (time5 - time4);
            log_debug("Finished particle velocity exchange between PEs (real time: %.4fs).\n", time5 - time4);
        }

        /*---------------------------------------------------------------
         * stress update ------------------------------------------------
         *---------------------------------------------------------------*/
        if (gv->FDORDER_TIME == 2) {

            switch (gv->WEQ) {
              case AC_ISO:     /* acoustic */
                  log_fatal("not yet implemented\n");
                  break;
              case AC_VTI:     /* acoustic VTI */
                  log_fatal("not yet implemented\n");
                  break;
              case AC_TTI:     /* acoustic TTI */
                  log_fatal("not yet implemented\n");
                  break;
              case EL_ISO:     /* elastic */
                  update_s_elastic_interior(nt, hc, mpm, mpw, gv);

                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_elastic_PML(nt, mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_elastic_abs(nt, mpm, mpw, gv);
                  }
                  break;
              case VEL_ISO:    /* viscoelastic */
                  update_s_visc_interior(nt, mpm, mpw, gv);
                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_visc_PML(nt, mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_visc_abs(mpm, mpw, gv);
                  }
                  break;
              case EL_VTI:     /* elastic VTI */
                  update_s_elastic_VTI_interior(nt, mpm, mpw, gv);
                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_elastic_VTI_PML(mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_elastic_vti_abs(mpm, mpw, gv);
                  }
#ifdef EBUG
                  debug_check_matrix(mpw->psxx, nt, gv->NX, gv->NY, 555, 0, "psxx");
                  debug_check_matrix(mpw->psyy, nt, gv->NX, gv->NY, 555, 0, "psyy");
                  debug_check_matrix(mpw->psxy, nt, gv->NX, gv->NY, 555, 0, "psxy");
#endif
                  break;
              case VEL_VTI:    /* viscoelastic VTI */
                  update_s_visc_VTI_interior(nt, mpm, mpw, gv);
                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_visc_VTI_PML(mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_visc_vti_abs(mpm, mpw, gv);
                  }
#ifdef EBUG
                  debug_check_matrix(mpw->psxx, nt, gv->NX, gv->NY, 666, 0, "psxx");
                  debug_check_matrix(mpw->psyy, nt, gv->NX, gv->NY, 666, 0, "psyy");
                  debug_check_matrix(mpw->psxy, nt, gv->NX, gv->NY, 666, 0, "psxy");
#endif
                  break;
              case EL_TTI:     /* elastic TTI */
                  update_s_elastic_TTI_interior(nt, mpm, mpw, gv);

                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_elastic_TTI_PML(mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_elastic_tti_abs(mpm, mpw, gv);
                  }
#ifdef EBUG
                  debug_check_matrix(mpw->psxx, nt, gv->NX, gv->NY, 777, 0, "psxx");
                  debug_check_matrix(mpw->psyy, nt, gv->NX, gv->NY, 777, 0, "psyy");
                  debug_check_matrix(mpw->psxy, nt, gv->NX, gv->NY, 777, 0, "psxy");
#endif
                  break;
              case VEL_TTI:    /* viscoelastic TTI */
                  update_s_visc_TTI_interior(nt, mpm, mpw, gv);
                  if (gv->FW) {
                      if (gv->ABS_TYPE == 1)
                          update_s_visc_TTI_PML(mpm, mpw, gv);
                      if (gv->ABS_TYPE == 2)
                          update_s_visc_tti_abs(mpm, mpw, gv);
                  }
#ifdef EBUG
                  debug_check_matrix(mpw->psxx, nt, gv->NX, gv->NY, 888, 0, "psxx");
                  debug_check_matrix(mpw->psyy, nt, gv->NX, gv->NY, 888, 0, "psyy");
                  debug_check_matrix(mpw->psxy, nt, gv->NX, gv->NY, 888, 0, "psxy");
#endif
                  break;
              case VAC_ISO:    /* viscoacoustic */
                  log_fatal("not yet implemented\n");
                  break;
              case VAC_VTI:    /* viscoacoustic VTI */
                  log_fatal("not yet implemented\n");
                  break;
              case VAC_TTI:    /* viscoacoustic TTI */
                  log_fatal("not yet implemented\n");
                  break;
              default:
                  log_fatal("Unknown WEQ.\n");
            }
        }

        if (gv->FDORDER_TIME == 4) {
            if (gv->L) {        /* viscoelastic */
                /* Not supported right now */
                update_s_visc_interior_4(nt, hc, mpm, mpw, gv);
                if (gv->FW) {
                    if (gv->ABS_TYPE == 1) {
                        update_s_visc_PML_4(nt, mpm, mpw, gv);
                    }
                    if (gv->ABS_TYPE != 1) {
                        update_s_visc_abs_4(nt, mpm, mpw, gv);
                    }
                }
                /* Shift memory variables one time step back */
                shift_var3(&(mpw->pp), &(mpw->pp_2), &(mpw->pp_3), &(mpw->pp_4));
                shift_var3(&(mpw->pr), &(mpw->pr_2), &(mpw->pr_3), &(mpw->pr_4));
                shift_var3(&(mpw->pq), &(mpw->pq_2), &(mpw->pq_3), &(mpw->pq_4));
            } else {            /* elastic */
                update_s_elastic_interior_4(nt, hc, mpm, mpw, gv);

                if (gv->FW) {
                    if (gv->ABS_TYPE == 1)
                        update_s_elastic_PML_4(nt, mpm, mpw, gv);
                    if (gv->ABS_TYPE != 1)
                        update_s_elastic_abs_4(nt, mpm, mpw, gv);
                }
            }
            /* Shift spatial derivatives from the velocity one time step back */
            shift_var2(&(mpw->vxx_1), &(mpw->vxx_2), &(mpw->vxx_3), &(mpw->vxx_4));
            shift_var2(&(mpw->vyy_1), &(mpw->vyy_2), &(mpw->vyy_3), &(mpw->vyy_4));
            shift_var2(&(mpw->vxy_1), &(mpw->vxy_2), &(mpw->vxy_3), &(mpw->vxy_4));
            shift_var2(&(mpw->vyx_1), &(mpw->vyx_2), &(mpw->vyx_3), &(mpw->vyx_4));
        }

        /* explosive source */
        if (gv->SOURCE_TYPE == 1)
            psource(nt, acq, mpw, gv);

        /* Applying free surface condition */
        if ((gv->FREE_SURF) && (gv->POS[2] == 0)) {
            if (gv->L)          /* viscoelastic */
                surface(1, hc, mpm, mpw, gv);
            else
                /* elastic */
                surface_elastic(1, hc, mpm, mpw, gv);
        }

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            time6 = MPI_Wtime();
            perf->time_av_s_update += (time6 - time5);
            log_debug("Starting stress exchange between PEs...\n");
        }

        /*---------------------------------------------------------------
         * -------- stress exchange between PEs --------
         *---------------------------------------------------------------*/
        exchange_s(mpw, gv);

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            time7 = MPI_Wtime();
            perf->time_av_s_exchange += (time7 - time6);
            log_debug("Finished stress exchange between PEs (real time: %.4fs).\n", time7 - time6);
        }

        /* store amplitudes at receivers in section-arrays */
        if ((gv->SEISMO) && (nt == lsamp) && (nt < gv->NT)) {

            seismo_ssg(lsamp, acq->recpos_loc, hc, mpm, mpw, gv);
            lsamp += gv->NDT;
        }

        /* WRITE SNAPSHOTS TO DISK */
        if ((gv->SNAP) && (nt == lsnap) && (nt <= gv->TSNAP2 / gv->DT)) {

            snap(nt, ++nsnap, hc, mpm, mpw, gv);
            lsnap = lsnap + iround(gv->TSNAPINC / gv->DT);
        }

        if ((gv->MPID == 0) && ((nt - 1) % gv->OUTNTIMESTEPINFO == 0)) {
            ++perf->infocounter;
            time8 = MPI_Wtime();
            perf->time_av_timestep += (time8 - time3);

            // when we reach this point, we have completed nt out of gv->NT time steps; however, the
            // time_av_timestep variable has only been updated every gv->OUTNTIMESTEPINFO time step;
            // use infocounter to calculate correct average
            log_info("Total real time for time step %d: %.4fs. Shot %d, time left: %.2lfs.\n", nt, time8 - time3,
                     ishot, (gv->NT - nt) * perf->time_av_timestep / (double)perf->infocounter);
        }
    }

    log_infoc(0, "Finished time stepping.\n");
}
