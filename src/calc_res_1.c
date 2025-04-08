
/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
*  calculate STF or residuals, observed data are read from file (+low-pass)
*  ----------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "fd.h"
#include "util.h"
#include "logging.h"

void calc_res_1(int ishot, int shot_id, st_seismogram *section, st_seismogram *section_obs, st_signals *signals,
              int ntr_loc, const float *finv, int nf, double *L2, int res_switch, int groupnum,
              Acq *acq, GlobVar *gv)
{

    switch (gv->ADJOINT_TYPE) {
      case 1:                  /* vx only */
          /* read seismic "observed" data (x-comp.) */
          readseis(shot_id, section_obs->vx, acq, ntr_loc, gv->NS, 1, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vx, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          if (gv->STFI_CALC) {
              if (gv->USE_TW) {
                  time_window(section->vx, section_obs->vx, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 1, gv);
              }
              /* STFI - Source time function inversion */
              stfi_calc(ishot, shot_id, section, section_obs, signals, ntr_loc, 1024, gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vx, section_obs->vx, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vx, section->vx, signals->sectionvxdiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          break;
      case 2:                  /* vy only */
          /* read seismic "observed" data (y-comp.) */
          readseis(shot_id, section_obs->vy, acq, ntr_loc, gv->NS, 2, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vy, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          if (gv->STFI_CALC) {
              if (gv->USE_TW) {
                  time_window(section->vy, section_obs->vy, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 1, gv);
              }
              /* STFI - Source time function inversion */
              stfi_calc(ishot, shot_id, section, section_obs, signals, ntr_loc, 1024, gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vy, section_obs->vy, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vy, section->vy, signals->sectionvydiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          break;
      case 3:                  /* vz only */
          /* read seismic "observed" data (z-comp.) */
          readseis(shot_id, section_obs->vz, acq, ntr_loc, gv->NS, 3, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vz, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          if (gv->STFI_CALC) {
              if (gv->USE_TW) {
                  time_window(section->vz, section_obs->vz, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 1, gv);
              }
              /* STFI - Source time function inversion */
              stfi_calc(ishot, shot_id, section, section_obs, signals, ntr_loc, 1024, gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vz, section_obs->vz, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vz, section->vz, signals->sectionvzdiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          break;
      case 4:                  /* particle velocities (vx + vy + vz) */
          /* read seismic "observed" data (x-comp.) */
          readseis(shot_id, section_obs->vx, acq, ntr_loc, gv->NS, 1, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vx, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vx, section_obs->vx, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vx, section->vx, signals->sectionvxdiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          /* read seismic "observed" data (y-comp.) */
          readseis(shot_id, section_obs->vy, acq, ntr_loc, gv->NS, 2, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vy, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          if (gv->STFI_CALC) {
              /* STFI - Source time function inversion */
              if (gv->USE_TW) {
                  time_window(section->vy, section_obs->vy, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 1, gv);
              }
              stfi_calc(ishot, shot_id, section, section_obs, signals, ntr_loc, 1024, gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vy, section_obs->vy, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vy, section->vy, signals->sectionvydiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          /* read seismic "observed" data (z-comp.) */
          readseis(shot_id, section_obs->vz, acq, ntr_loc, gv->NS, 3, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->vz, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->vz, section_obs->vz, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->vz, section->vz, signals->sectionvzdiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          break;
      case 5:                  /* pressure only */
          /* read seismic "observed" data (pressure) */
          readseis(shot_id, section_obs->p, acq, ntr_loc, gv->NS, 4, gv);
          if (gv->FILT == 1) {
              filt_seis(section_obs->p, ntr_loc, gv->NS, finv[nf - 1], gv);
          }
          if (gv->STFI_CALC) {
              /* STFI - Source time function inversion */
              if (gv->USE_TW) {
                 time_window(section->p, section_obs->p, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 1, gv);
              }
              stfi_calc(ishot, shot_id, section, section_obs, signals, ntr_loc, 1024, gv);
          }
          /* calculate residuals and L2 Norm */
          if (res_switch) {
              if (gv->USE_TW) {
                  time_window(section->p, section_obs->p, acq->recpos_loc, gv->NTRG, ntr_loc, shot_id, groupnum, 0, gv);
              }
              residual(section_obs->p, section->p, signals->sectionpdiff, ntr_loc, gv->NS, L2, finv, nf, gv);
          }
          break;
    }
    gv->STFI_CALC = 0;
}
