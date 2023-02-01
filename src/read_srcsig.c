
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
 *   Read external source wavelet (ASCII or SU format)                        
 *------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "read_srcsig.h"
#include "read_su.h"
#include "macros.h"
#include "globvar_struct.h"
#include "logging.h"
#include "su_struct.h"

float *read_srcsig(int *ns, const GlobVar * gv)
{
    const char *sigfile = gv->SIGNAL_FILE;
    const char *sig_filetype[] = { "unknown", "ASCII", "SU" };

    char buffer[STRING_SIZE], *p = NULL, *q = NULL;
    int filetype = 0, i = 0, j = 0;
    float *sig = NULL, dt = 0.f;
    FILE *filep = NULL;
    bool dt_found = false;
    size_t nt = 0;
    short sudelrt = 0;
    unsigned short suns = 0, sudt = 0;
    SUhead suhead;

    log_infoc(0, "Source signature file: %s\n", gv->SIGNAL_FILE);

    // try to automatically determine what type of file it is based on file suffix
    if (STRSTRCOMP(sigfile, ".su")) {
        filetype = 2;
    } else if (STRSTRCOMP(sigfile, ".dat") || STRSTRCOMP(sigfile, ".txt") || STRSTRCOMP(sigfile, ".asc")) {
        filetype = 1;
    } else {
        filetype = 0;
    }

    log_debugc(0, "Automatically determined file type seems to be '%s'.\n", sig_filetype[filetype]);

    *ns = 0;

    if (filetype < 2) {

                    /***** try ASCII *****/

        filep = fopen(sigfile, "r");
        if (!filep)
            log_fatal("Could not open file %s for reading.\n", sigfile);

        i = 0;
        while (fgets(buffer, STRING_SIZE, filep)) {
            ++i;
            log_debugc(0, "Read line %d, content: %s", i, buffer);
            p = buffer;
            while (isspace((unsigned char)*p))
                ++p;
            if (*p == '\0' || (*p == '\n') || (*p == '\r'))
                continue;
            if (*p == '#') {
                ++p;
                if (*p == '\0' || (*p == '\n') || (*p == '\r'))
                    continue;
                q = STRSTRCOMP(p, "DT");
                if (!q)
                    continue;
                q = STRSTRCOMP(p, "=");
                if (!q)
                    continue;
                ++q;
                if (*q == '\0' || (*p == '\n') || (*p == '\r'))
                    continue;
                while (isspace((unsigned char)*q))
                    ++q;
                j = sscanf(q, "%e", &dt);
                if (1 == j) {
                    log_debugc(0, "Source signature sampling interval: %fs\n", dt);
                    dt_found = true;
                }
                continue;
            }
            ++(*ns);
        }

        log_infoc(0, "Number of samples in source signature file: %d\n", *ns);
        if (0 == *ns && 0 == gv->MPID)
            log_fatal("Source signature file did not contain valid samples.\n");

        sig = (float *)calloc(sizeof(float), *ns);
        if (!sig)
            log_fatal("Could not allocate source signature buffer.\n");

        rewind(filep);

        i = 0;
        while (fgets(buffer, STRING_SIZE, filep)) {
            p = buffer;
            while (isspace((unsigned char)*p))
                ++p;
            if (*p == '\0' || (*p == '\n') || (*p == '\r') || (*p == '#'))
                continue;
            j = sscanf(p, "%e", &(sig[i]));
            if (0 == j)
                log_warnc(0, "Scan error while reading source signature. Check result!\n");
            ++i;
        }

    } else if (2 == filetype) {

                            /***** try SU *****/

        filep = fopen(sigfile, "rb");
        if (!filep)
            log_fatal("Could not open SU file %s for reading.\n", sigfile);

        nt = su_get_nt(filep, &suns, &sudt, &sudelrt);
        if (nt > 1)
            log_warnc(0, "SU file %s contains more than one trace. Reading only first trace.\n", sigfile);

        *ns = suns;
        log_debugc(0, "SU file %s contains %ld trace(s).\n", sigfile, nt);
        log_debugc(0, "Samples per trace: %d, dt=%d, delrt=%d.\n", *ns, (int)sudt, (int)sudelrt);

        sig = (float *)calloc(sizeof(float), *ns);
        if (!sig)
            log_fatal("Could not allocate source signature buffer.\n");

        su_read_trace(filep, 0, suns, true, &suhead, sig);

        if (sudt > 0) {
            dt = sudt * 1e-6;
            dt_found = true;
        }
        if (sudelrt != 0) {
            if (0 == gv->MPID)
                log_fatal("Delay recording time != 0 not supported for source signature.\n");
        }

    } else {
        if (0 == gv->MPID)
            log_fatal("Unknown filetype in source signature read function.\n");
    }

    if (filep)
        fclose(filep);

    if (dt_found) {
        if (dt != gv->DT && 0 == gv->MPID) {
            log_error("dt source signature: %fs, FD step size: %fs\n", dt, gv->DT);
            log_fatal("Sampling interval of source signature does not match FD step size!\n");
        } else {
            log_infoc(0, "Sampling interval of source signature matches FD step size.\n");
        }
    } else {
        log_warnc(0, "No sampling interval available for source signature. Make sure dt is correct.\n");
    }

    return sig;
}

#ifdef READ_SRCSIG_MAIN

// test main program //
int main(int argc, char **argv)
{
    log_init(NULL);
    log_set_level(LOG_DEBUG);

    if (argc < 2) {
        log_fatal("No source wavelet input file name given.\n");
    }

    GlobVar gv;
    strncpy(&(gv.SIGNAL_FILE[0]), argv[1], 255);
    gv.DT = 0.0001;
    log_info("Input file: %s\n", gv.SIGNAL_FILE);

    int nts;
    float *srcsig = read_srcsig(&nts, &gv);

    for (int i = 0; i < nts; ++i) {
        log_std("%f\n", srcsig[i]);
    }

    free(srcsig);

    log_finalize();

    exit(EXIT_SUCCESS);
}

#endif
