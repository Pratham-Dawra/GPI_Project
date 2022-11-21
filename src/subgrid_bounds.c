/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/* $Id: subgrid_bounds.c 819 2015-04-17 11:07:06Z tmetz $ */
/*---------------------------------------------------------------
calculate for-loop bounds for each process/subgrid

   [1] [2]                 [3] [4] gx
 [1]|---|-------------------|---|
    |   |                   |   |
 [2]|---|-------------------|---|
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
    |   |                   |   |
 [3]|---|-------------------|---|
    |   |                   |   |
 [4]|---|-------------------|---|

 gy

Width of Absorbing Boundary is gv->FW Gridpoints.
nx1,nx2,ny1,ny2 are first are the bounds for each subgrid without
absorbing boundaries.
In case of number of processes NPROCX=1 and NPROCY=1 all boundaries are in one
grid.
For NPROCX>2 and NPROCY>2 there are 9 different subgrids.
gx[2],gy[2] is the last gridpoint (FW) of the left boundary 
gx[3],gy[3] is the last gridpoint of the midpart  gx[3]=(nx2-FW),gy[3]=(ny2-FW)

Free Surface and Periodic Boundary conditions alter the bounds (see end of this file)
*--------------------------------------------------------------*/
#include "fd.h"


void subgrid_bounds ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, GlobVar *gv )
//void subgrid_bounds ( int nx1, int nx2, int ny1, int ny2)
{

	/* GRID */
	switch ( gv->NPROCY ) {

	case 1: /*case gv->NPROCY=1 */

		switch ( gv->NPROCX ) {
		case 1:     /* Case gv->NPROCX=1 */
			gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
			gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			
			break;

		case 2:   /* Case gv->NPROCX=2 */
			if ( gv->POS[1]==0 ) { /* left */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( gv->POS[1]== gv->NPROCX-1 ) { /* right */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			break;

		default: /* Case gv->NPROCX>2 */

			if ( gv->POS[1]==0 ) { /* left */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}

			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!= gv->NPROCX-1 ) ) { /* mid */

				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( gv->POS[1]== gv->NPROCX-1 ) { /* right */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			break;
		}
		break;

	case 2: /*case gv->NPROCY=2 */
		switch ( gv->NPROCX ) {
		case 1:   /* Case gv->NPROCX=1 */
			if ( gv->POS[2]==0 ) { /* top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( gv->POS[2]==gv->NPROCY-1 ) { /* bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			break;

		case 2: /* Case gv->NPROCX=2 */
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]==0 ) ) { /* corner left-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]== gv->NPROCY-1 ) ) { /* corner left-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /* corner right-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* corner right-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			break;

		default: /* Case gv->NPROCX>2 */

			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]==0 ) ) { /* corner left-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]== gv->NPROCY-1 ) ) { /* corner left-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /* corner right-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* corner right-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!=gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /*top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!=gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}


			break;
		} /* end of switch gv->NPROCX */
		break;

	default: /*case gv->NPROCY>2 */
		switch ( gv->NPROCX ) {
		case 1:   /* Case gv->NPROCX=1 */
			if ( gv->POS[2]==0 ) { /* top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			if ( ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* mid */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			if ( gv->POS[2]==gv->NPROCY-1 ) { /* bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			break;

		case 2: /* Case gv->NPROCX=2 */
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]==0 ) ) { /* corner left-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]== gv->NPROCY-1 ) ) { /* corner left-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /* corner right-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* corner right-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* left */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* right */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			break;

		default: /* Case gv->NPROCX>2 */

			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]==0 ) ) { /* corner left-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]== gv->NPROCY-1 ) ) { /* corner left-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /* corner right-top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* corner right-bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			if ( ( gv->POS[1]==0 ) && ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* left */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=gv->FW, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]==gv->NPROCX-1 ) && ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* right */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2-gv->FW, gx[4]=nx2;
			}

			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!=gv->NPROCX-1 ) && ( gv->POS[2]==0 ) ) { /*top */
				gy[1]=ny1, gy[2]=gv->FW, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}
			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!=gv->NPROCX-1 ) && ( gv->POS[2]==gv->NPROCY-1 ) ) { /* bottom */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2-gv->FW, gy[4]=ny2;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}

			if ( ( gv->POS[1]!=0 ) && ( gv->POS[1]!=gv->NPROCX-1 ) && ( gv->POS[2]!=0 ) && ( gv->POS[2]!=gv->NPROCY-1 ) ) { /* center */
				gy[1]=ny1, gy[2]=ny1-1, gy[3]=ny2, gy[4]=ny2-1;
				gx[1]=nx1, gx[2]=nx1-1, gx[3]=nx2, gx[4]=nx2-1;
			}
			break;
		} /*end of switch (gv->NPROCX) */
		break;
	} /*end of switch (gv->NPROCY) */

	if ( gv->FREE_SURF==1 )
		if ( gv->POS[2] == 0 )
			gy[2]=ny1-1;

		
        /* Periodic Boundary condition */
	if ( gv->BOUNDARY==1 ) { 
		if ( gv->POS[1] == 0 )
			gx[2]=nx1-1;
		if ( gv->POS[1] == gv->NPROCX-1 )
			gx[3]=nx2,gx[4]=nx2-1;
	}
}



