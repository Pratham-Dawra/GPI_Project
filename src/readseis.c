
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
 * reads observed seismograms from files;
 * this function is not thread-safe
 *------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "fd.h"
#include "macros.h"
#include "logging.h"
#include "read_su.h"
#include "su_struct.h"
#include "interpol.h"

static bool b_init = false;
static float *tout = NULL;

void readseis(int ishot, float **section, const AcqVar *acq, int ntr, int ns, int comp, const GlobVar *gv)
{

    char data[2*STRING_SIZE], shotnum[16];
    unsigned short suns, sudt;
    short sudelrt;

    const char *component[] = { "vx", "vy", "vz" , "p"};

    if (comp > 4) {
	log_fatal("Unknown component %d requested in readseis.\n", comp);
    }

    if (gv->RUN_MULTIPLE_SHOTS) {
	sprintf(shotnum, ".shot%d", ishot);
    } else {
	sprintf(shotnum, "%s", "");
    }

    sprintf(data, "%s_%s.su%s", gv->SEIS_OBS_FILE, component[comp-1], shotnum);
	if(gv->MPID==0) log_info("Reading file %s.\n", data);

    if (ntr > 0) {
        bool b_interp = false;
		FILE *fpdata = fopen(data, "r");
        if (!fpdata) {
            log_fatal("Could not open file %s for reading.", data);
        }
        size_t sunt  = su_get_nt(fpdata, &suns, &sudt, &sudelrt);
        if (sunt < (unsigned)ntr) {
	    log_fatal("File %s: fewer traces (%d) in the file than we need to read (%d).\n", 
		      data, sunt, ntr);
        }
	
	// use dt from SU file unless user gave us an explicit dt in JSON file
	float dt = sudt*1e-6;
	if (gv->SEIS_OBS_DT > 0.0) {
	    dt = gv->SEIS_OBS_DT;
	}
	if (dt != gv->DT) {
	  // interpolation required if time step interval and seismogram sampling interval differ
	  b_interp = true;
	}
	
	int k = 0;
	if (b_interp) {
            if (!b_init) {
	        tout = (float*)malloc(gv->NS*sizeof(float));
		for (int i=0; i<gv->NS; i++) {
	            tout[i] = i*gv->DT;
	        }
		b_init = true;
	    }
	    float *trace = (float*)malloc(suns*sizeof(float));
	    // SUhead header;
	    for (size_t i = 1; i <= sunt; i++) {
	        if (1 == acq->recswitch[i]) {
	          // su_read_trace(fpdata, i - 1, (unsigned)ns, true, &header, &section[++k][0]);
	          su_read_trace(fpdata, i - 1, suns, true, NULL, trace);
		  sinc_ipol(suns, dt, 0.0f, trace, gv->NS, tout, &section[++k][0]);
	        }
	    }   
	    free(trace);
	    // we deliberately leak tout at this point to avoid having to re-create
	    // the time vector over and over again; as gv->NS does not change during
	    // a run, it could be avoided by making the vector part of the st_seismogram
	    // structure
	} else {
	    if (suns < (unsigned)ns) {
                log_fatal("File %s: fewer number of samples (%d) in the file than we need to read (%d).\n",
		          data, suns, ns);
	    } else if (suns > (unsigned)ns) {
	        log_warn("File %s: more samples (%d) in the file than required (%d).\n", data, suns, ns);
	    }
	    // SUhead header;
	    for (size_t i = 1; i <= sunt; i++) {
	        if (1 == acq->recswitch[i]) {
	          // su_read_trace(fpdata, i - 1, (unsigned)ns, true, &header, &section[++k][0]);
	          su_read_trace(fpdata, i - 1, (unsigned)ns, true, NULL, &section[++k][0]);
	        }
	    }
	}
	if (k != ntr) log_fatal("k: %d, ntr: %d", k, ntr);
	assert(k == ntr); // cross-check: must have read ntr traces
        fclose(fpdata);
    }
}
