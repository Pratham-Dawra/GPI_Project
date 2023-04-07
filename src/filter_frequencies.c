
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

/* 
   Reading filter frequencies for the time domain filter from FREQ_FILE.
   
*/

#include "fd.h"
#include "logging.h"

void filter_frequencies(GlobVarInv *vinv)
{
    FILE *freqf;
    char buffer[STRING_SIZE];
    float buffer_f;
    int l;

    /*---------------------------- open FREQ_FILE ---------------------------*/
    if ((freqf = fopen(vinv->FREQ_FILE, "r")) == NULL)
        log_fatal("Freqency file could not be opened!\n");

    log_info("Reading frequencies from file: %s\n", vinv->FREQ_FILE);

    /* Count how many lines the work flow file has */
    while (fgets(buffer, STRING_SIZE, freqf)) {
        /* checks if the line contains a '#'character which indicates a comment line,
         * and if the line contains a float */
        if ((sscanf(buffer, "%f", &buffer_f) == 1) && buffer_f != 0) {
            ++(vinv->NFREQ);
        }
    }

    /*-------------------alocate frequency array-----------------------------*/
    vinv->F_LOW_PASS = vector(1, vinv->NFREQ);

    /*------------------Read frequencies from FREQ_FILE----------------------*/
    rewind(freqf);
    l = 1;
    while (fgets(buffer, STRING_SIZE, freqf)) {
        //sscanf(buffer, "%s", bufferstring);
        if ((sscanf(buffer, "%f", &(vinv->F_LOW_PASS)[l]) == 1) && (vinv->F_LOW_PASS[l] != 0)) {
            if (l > vinv->NFREQ) {
                log_fatal
                    ("read_workflow.c: buffer not large enough to store all workflow parameters - programming error.\n");
            }
            ++l;
        }
    }

    fclose(freqf);

    return;
}
