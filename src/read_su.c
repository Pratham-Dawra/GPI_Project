
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
 *   Read data in (native-endian) SU format
 *------------------------------------------------------------------------*/

#include "read_su.h"
#include "su_gather.h"
#include "su_struct.h"
#include "logging.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* determine SU file size */
off_t su_filesize(FILE * filep)
{
    int ierr = fseek(filep, 0, SEEK_END);
    if (ierr != 0)
        log_fatal("Could not seek to end of SU file.\n");
    return ftell(filep);
}

size_t su_get_nt(FILE * filep, unsigned short *ns, unsigned short *dt, short *delrt)
{
    int ierr = 0;
    SUhead header;
    size_t tracelen, nt;
    off_t filesize = su_filesize(filep);

    *ns = 0;
    *dt = 0;
    *delrt = 0;

    log_debug("Length of SU file in bytes: %ld, size of SU header in bytes: %d\n", filesize, sizeof(SUhead));

    ierr = fseek(filep, 0, SEEK_SET);
    if (ierr != 0)
        log_fatal("Could not seek to beginning of SU file.\n");

    ierr = fread(&header, sizeof(SUhead), 1, filep);
    if (ierr != 1)
        log_fatal("Short read in SU file.\n");

    if (header.ns < 1)
        log_fatal("SU file has no number of samples (ns) in first SU header.\n");

    *ns = header.ns;
    *dt = header.dt;
    *delrt = header.delrt;

    tracelen = sizeof(SUhead) + header.ns * sizeof(float);

    if (fmod((double)filesize, (double)tracelen)) {
        log_warn("Lenght mismatch in SU file. Filesize: %ld bytes, trace length: %ld bytes. Check SU file.\n", filesize,
                 tracelen);
    }

    nt = filesize / tracelen;

    ierr = fseek(filep, 0, SEEK_SET);
    if (ierr != 0)
        log_fatal("Could not seek to beginning of SU file.\n");

    return nt;
}

int su_read_trace(FILE * filep, size_t n, unsigned short ns, bool b_seek, SUhead * header, float *data)
{
    int ierr = 0;
    size_t tracelen = 0;
    off_t fileoffset = 0;

    // seek to correct trace position if requested
    if (b_seek) {
        tracelen = sizeof(SUhead) + ns * sizeof(float);
        fileoffset = n * tracelen;

        log_debug("File offset to read trace %ld: %ld bytes\n", n, fileoffset);

        ierr = fseek(filep, fileoffset, SEEK_SET);
        if (ierr != 0)
            log_fatal("Could not seek to offset %ld of SU file.\n", fileoffset);
    }
    // read SU trace header if requested
    if (header) {
        ierr = fread(header, sizeof(SUhead), 1, filep);
        if (ierr != 1)
            log_fatal("Short read in SU file while reading header.\n");
    } else {
        ierr = fseek(filep, sizeof(SUhead), SEEK_CUR);
        if (ierr != 0)
            log_fatal("Could not skip over trace header in SU file.\n");
    }

    // read SU trace data
    ierr = fread(data, sizeof(float), ns, filep);
    if (ierr != ns)
        log_fatal("Short read in SU file while reading data.\n");

    return 1;
}

size_t su_read_file(FILE * filep, SUgather * gather)
{
    size_t i, total = 0;
    unsigned short ns, dt;
    short delrt;

    size_t nt = su_get_nt(filep, &ns, &dt, &delrt);

    malloc_SUgather(gather, nt, ns);

    for (i = 0; i < nt; ++i) {
        total += su_read_trace(filep, i, ns, false, &(gather->header[i]), &(gather->data[i][0]));
    }

    return total;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef READ_SU_MAIN

// test main program //
int main(int argc, char **argv)
{
    FILE *filep = NULL;

    log_init(NULL);
    log_set_level(LOG_DEBUG);

    if (argc < 2) {
        log_fatal("No SU input file name given.\n");
    }

    const char *sufile = argv[1];
    log_info("Input file: %s\n", sufile);

    filep = fopen(sufile, "rb");
    if (!filep)
        log_fatal("Could not open SU file %s for reading.\n", sufile);

    SUgather gather;
    size_t nt_read = su_read_file(filep, &gather);

    if (nt_read != gather.nt)
        log_error("Internal 'number of traces' mismatch. Programming error.\n");

    log_info("Read gather with %ld traces and %d samples per trace.\n", gather.nt, (int)gather.ns);

    free_SUgather(&gather);

    if (filep)
        fclose(filep);

    log_finalize();

    exit(EXIT_SUCCESS);
}

#endif

///////////////////////////////////////////////////////////////////////////////
