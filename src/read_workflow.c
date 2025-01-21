
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
 *  Read extern workflow                          
 *
 *
 *  If you want to adjust the workflow, I think it is not necessary to
 *  modify this file. Have a look at apply_workflow.c and adjust
 *  WORKFLOW_MAX_VAR in fd.h.
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include "macros.h"

void read_workflow(GlobVar *gv, GlobVarInv *vinv)
{
    /* workflow is a pointer to a pointer, keep care... */

    /* intern variables */
    int l;
    char buffer[STRING_SIZE], bufferstring[10];
    FILE *fwork;
    /* WORKFLOW_MAX_VAR is set in macros.h */

    if (gv->MPID == 0)
        log_info("Reading workflow from file: %s\n", vinv->FILE_WORKFLOW);

    /* Open Workflow file */
    if ((fwork = fopen(vinv->FILE_WORKFLOW, "r")) == NULL) {
        log_fatal("Workflow file  %s could not be opened.\n", vinv->FILE_WORKFLOW);
    }

    /* Count how many lines the work flow file has */
    while (fgets(buffer, STRING_SIZE, fwork)) {
        sscanf(buffer, "%s", bufferstring);
        /* checks if the line contains a '#'character which indicates a comment line,
         * and if the reading of a string was successful, which is not the case for an empty line */
        if ((strchr(buffer, '#') == 0) && (sscanf(buffer, "%s", bufferstring) == 1)) {
            ++(vinv->WORKFLOW_LINES);
        }
    }

    if (vinv->WORKFLOW_LINES == 0) {
        if (gv->MPID == 0)
            log_warn("Could not determine number of workflow parameter sets; assuming %d.\n", vinv->WORKFLOW_LINES);
    } else {
        vinv->WORKFLOW_LINES -= 1;  // Don't count header line.
    }

    /* Allocate memory for workflow */
    vinv->WORKFLOW = matrix(1, vinv->WORKFLOW_LINES, 1, WORKFLOW_MAX_VAR);
    if ((vinv->WORKFLOW) == NULL)
        log_fatal("Was not able to allocate memory for WORKFLOW buffer!");

    rewind(fwork);
    fgets(vinv->WORKFLOW_HEADER, STRING_SIZE, fwork);   /* Read header from first line */
    l = 0;
    while (fgets(buffer, STRING_SIZE, fwork)) {
        sscanf(buffer, "%s", bufferstring);
        if ((strchr(buffer, '#') == 0) && (sscanf(buffer, "%s", bufferstring) == 1)) {
            ++l;
            if (l > vinv->WORKFLOW_LINES) {
                log_fatal
                    ("read_workflow.c: buffer not large enough to store all workflow parameters - programming error.\n");
            }
            sscanf(buffer, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &(vinv->WORKFLOW)[l][1], &(vinv->WORKFLOW)[l][2],
                   &(vinv->WORKFLOW)[l][3], &(vinv->WORKFLOW)[l][4], &(vinv->WORKFLOW)[l][5], &(vinv->WORKFLOW)[l][6],
                   &(vinv->WORKFLOW)[l][7], &(vinv->WORKFLOW)[l][8], &(vinv->WORKFLOW)[l][9], &(vinv->WORKFLOW)[l][10],
                   &(vinv->WORKFLOW)[l][11], &(vinv->WORKFLOW)[l][12], &(vinv->WORKFLOW)[l][13],
                   &(vinv->WORKFLOW)[l][14]);
            /* set EPRECOND to the maximum in workflow */
            if (vinv->EPRECOND_MAX < (vinv->WORKFLOW)[l][8]) vinv->EPRECOND_MAX = (vinv->WORKFLOW)[l][8];
        }
    }
    fclose(fwork);
}
