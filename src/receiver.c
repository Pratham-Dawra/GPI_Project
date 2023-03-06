
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

/*----------------------------------------------------------------------
 * compute receiver positions or read them from external file
 *----------------------------------------------------------------------*/

#include "fd.h"
#include "logging.h"
#include <math.h>
#include <stdbool.h>

int **receiver(GlobVar *gv)
{
    int **recpos1 = NULL, **recpos = NULL;
    int itr = 1, itr1 = 0, itr2 = 0, recflag = 0, n, i, j;
    float nxrec = 0, nyrec = 0;
    float nxrec1, nxrec2, nyrec1, nyrec2;
    float xrec, yrec;
    bool testbuff1, testbuff2, testbuff3;
    char bufferstring[10], buffer[STRING_SIZE];
    FILE *fpr = NULL;

    if (gv->MPID == 0) {
        log_info("------------------------- Receiver positions ----------------\n");
        if (gv->READREC) {      /* read receiver positions from file */
            log_info("Reading receiver positions from file %s.\n", gv->REC_FILE);
            fpr = fopen(gv->REC_FILE, "r");
            if (fpr == NULL)
                log_fatal("Receiver file %s could not be opened.\n", gv->REC_FILE);
            gv->NTRG = 0;

            /* counts the number of receivers in the receiver file */
            while (fgets(buffer, STRING_SIZE, fpr)) {
                testbuff1 = strchr(buffer, '#');
                testbuff2 = strchr(buffer, '%');
                testbuff3 = sscanf(buffer, "%s", bufferstring) == 1;

                /* checks if the line contains a '%' or '#' character which indicates a
                 * comment line, and if the reading of a string was successful, 
                 * which is not the case for an empty line */
                if (((testbuff1 == 1 || testbuff2 == 1) == 0) && testbuff3 == 1)
                    ++(gv->NTRG);
            }

            rewind(fpr);

            recpos1 = imatrix(1, 3, 1, gv->NTRG);
            for (itr = 1; itr <= gv->NTRG; itr++) {
                fscanf(fpr, "%f%f\n", &xrec, &yrec);
                recpos1[1][itr] = (int)floor((xrec + gv->REFREC[1]) / gv->DH) + 1;
                recpos1[2][itr] = (int)floor((yrec + gv->REFREC[2]) / gv->DH) + 1;
                recpos1[3][itr] = itr;
            }
            fclose(fpr);

            /* check if more than one receiver is located at the same gridpoint */
            for (itr = 1; itr <= (gv->NTRG - 1); itr++)
                for (itr1 = itr + 1; itr1 <= gv->NTRG; itr1++)
                    if ((recpos1[1][itr] == recpos1[1][itr1]) && (recpos1[2][itr] == recpos1[2][itr1]))
                        recpos1[1][itr1] = -(++recflag);

            recpos = imatrix(1, 3, 1, gv->NTRG - recflag);
            for (itr = 1; itr <= gv->NTRG; itr++)
                if (recpos1[1][itr] > 0) {
                    recpos[1][++itr2] = recpos1[1][itr];
                    recpos[2][itr2] = recpos1[2][itr];
                    recpos[3][itr2] = recpos1[3][itr];
                }

            int ntr_orig = gv->NTRG;
            gv->NTRG = itr2;
            if (recflag > 0) {
                log_warn("Co-located receivers positions; reduced no. of receivers from %d to %d.\n",
                         ntr_orig, gv->NTRG);
            }

            free_imatrix(recpos1, 1, 3, 1, gv->NTRG);
        } else if (gv->REC_ARRAY > 0) {
            log_info("Generating receiver planes as specified in input file.\n");

            gv->NTRG = (1 + (gv->NXG - 2 * gv->FW) / gv->DRX) * gv->REC_ARRAY;
            recpos = imatrix(1, 3, 1, gv->NTRG);
            itr = 0;
            for (n = 0; n <= gv->REC_ARRAY - 1; n++) {
                j = iround((gv->REC_ARRAY_DEPTH + gv->REC_ARRAY_DIST * (float)n) / gv->DH);
                for (i = gv->FW; i <= gv->NXG - gv->FW; i += gv->DRX) {
                    itr++;
                    recpos[1][itr] = i;
                    recpos[2][itr] = j;
                    recpos[3][itr] = itr;
                }
            }
        } else {                /* straight horizontal or vertical line of receivers */
            log_info("Receiver positions specified by json parameter file.\n");
            nxrec1 = (int)floor(gv->XREC1 / gv->DH) + 1;    /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
            nyrec1 = (int)floor(gv->YREC1 / gv->DH) + 1;    /* the positions of the first and last receiver */
            nxrec2 = (int)floor(gv->XREC2 / gv->DH) + 1;    /* in gridpoints */
            nyrec2 = (int)floor(gv->YREC2 / gv->DH) + 1;

            /* only 1 receiver */
            if ((nxrec2 - nxrec1) == 0 && (nyrec2 - nyrec1) == 0) {
                log_info("A single receiver position determined from json file.\n");
                gv->NTRG = 1;
                recpos = imatrix(1, 3, 1, gv->NTRG);
                recpos[1][1] = nxrec1;
                recpos[2][1] = nyrec1;
                recpos[3][1] = 1;
            } else if ((nyrec2 - nyrec1) == 0) {
                log_info("A horizontal receiver line (in x-direction) determined from json file.\n");
                gv->NTRG = (int)floor((nxrec2 - nxrec1) / gv->NGEOPH) + 1;
                recpos = imatrix(1, 3, 1, gv->NTRG);
                nxrec = nxrec1;
                for (n = 1; n <= gv->NTRG; n++) {
                    nyrec = nyrec1 + ((nyrec2 - nyrec1) / (nxrec2 - nxrec1) * (nxrec - nxrec1));
                    itr = (int)floor((nxrec - nxrec1) / gv->NGEOPH) + 1;
                    recpos[1][itr] = nxrec;
                    recpos[2][itr] = nyrec;
                    recpos[3][itr] = itr;
                    nxrec += gv->NGEOPH;
                }
            } else if ((nxrec2 - nxrec1) == 0) {
                log_info("A vertical receiver line (in y-direction) determined from json file.\n");
                gv->NTRG = (int)floor((nyrec2 - nyrec1) / gv->NGEOPH) + 1;
                recpos = imatrix(1, 3, 1, gv->NTRG);
                nyrec = nyrec1;
                for (n = 1; n <= gv->NTRG; n++) {
                    nxrec = nxrec1 + ((nxrec2 - nxrec1) / (nyrec2 - nyrec1) * (nyrec - nyrec1));
                    itr = (int)floor((nyrec - nyrec1) / gv->NGEOPH) + 1;
                    recpos[1][itr] = nxrec;
                    recpos[2][itr] = nyrec;
                    recpos[3][itr] = itr;
                    nyrec += gv->NGEOPH;
                }
            } else {
                /* arbitrary geophone-line */
                log_error("No horizontal or vertical receiver line is specified in the json file.\n");
                log_error
                    ("In order to define an arbitrary receiver line, please make use of an external receiver file (READREC=1).\n");
                log_fatal("Error in specifying receiver coordinates in the json file!\n");
            }
        }                       /* end of if receivers specified in input file */
        log_info("Number of receiver positions: %d\n", gv->NTRG);

        for (itr = 1; itr <= gv->NTRG; itr++) {
            if (recpos[1][itr] == 0) recpos[1][itr] = 1;
            if (recpos[2][itr] == 0) recpos[2][itr] = 1;
        }
    }
    /* End of if(gv->MPID==0) */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(gv->NTRG), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (gv->MPID != 0)
        recpos = imatrix(1, 3, 1, gv->NTRG);
    MPI_Bcast(&recpos[1][1], (gv->NTRG) * 3, MPI_INT, 0, MPI_COMM_WORLD);

    if (gv->MPID == 0) {
        if (gv->NTRG > 50)
            log_warn("The following table is quite large (%d lines); only printing the first 50 entries!\n", gv->NTRG);
        log_info("Receiver positions in the global model (indices are one-based):\n");
        log_info("x (gridpts)  y (gridpts)      x (in m)      y (in m)\n");
        log_info("-----------  -----------  ------------  ------------\n");

        int maxprint = gv->NTRG;
        if (maxprint > 50)
            maxprint = 50;

        for (itr = 1; itr <= maxprint; itr++) {
            log_info("%11d  %11d  %12.2f  %12.2f\n",
                     recpos[1][itr], recpos[2][itr], (recpos[1][itr]-1) * gv->DH, (recpos[2][itr]-1) * gv->DH);
        }
    }

    return recpos;
}
