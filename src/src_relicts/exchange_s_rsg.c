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
 *   write values of dynamic field variables at the edges of the
 *   local grid into buffer arrays (which will be exchanged between
 *   processes)
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_s_rsg(float ** vx, float ** vy, float ** vz,
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, GlobVar *gv){

	int i, j;
	MPI_Status  status;



	
	if (gv->POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=gv->NX;i++){

			/* storage of top of local volume into buffer */
			buffertop_to_bot[i][1]  =  vx[1][i];
			buffertop_to_bot[i][2]  =  vy[1][i];
			buffertop_to_bot[i][3]  =  vz[1][i];
			
	}


	if (gv->POS[2]!=gv->NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=gv->NX;i++){
		
			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][1]  =  vx[gv->NY][i];
			bufferbot_to_top[i][2]  =  vy[gv->NY][i];
			bufferbot_to_top[i][3]  =  vz[gv->NY][i];

	}

	
  	 /* send and reveive values for points at inner boundaries */

	MPI_Sendrecv_replace(&buffertop_to_bot[1][1],gv->NX*3,MPI_FLOAT,gv->INDEX[3],gv->TAG5,gv->INDEX[4],gv->TAG5,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&bufferbot_to_top[1][1],gv->NX*3,MPI_FLOAT,gv->INDEX[4],gv->TAG6,gv->INDEX[3],gv->TAG6,MPI_COMM_WORLD,&status);


	
	if (gv->POS[2]!=gv->NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=gv->NX;i++){

			vx[gv->NY+1][i] = buffertop_to_bot[i][1];
			vy[gv->NY+1][i] = buffertop_to_bot[i][2];
			vz[gv->NY+1][i] = buffertop_to_bot[i][3];
			
	}

	if (gv->POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=gv->NX;i++){

			vx[0][i] = bufferbot_to_top[i][1];
			vy[0][i] = bufferbot_to_top[i][2];
			vz[0][i] = bufferbot_to_top[i][3];
			
	}
				
	
	
	/* exchange if periodic boundary condition is applied */
	if ((gv->BOUNDARY) || (gv->POS[1]!=0))	
	for (j=0;j<=gv->NY+1;j++){


			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][1] =  vx[j][1];
			bufferlef_to_rig[j][2] =  vy[j][1];
			bufferlef_to_rig[j][3] =  vz[j][1];


	}


														/* no exchange if periodic boundary condition is applied */
	if ((gv->BOUNDARY) || (gv->POS[1]!=gv->NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=0;j<=gv->NY+1;j++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][1] =  vx[j][gv->NX];
			bufferrig_to_lef[j][2] =  vy[j][gv->NX];
			bufferrig_to_lef[j][3] =  vz[j][gv->NX];


	}

  
  	 /* send and reveive values for points at inner boundaries */

 	MPI_Sendrecv_replace(&bufferlef_to_rig[0][1],(gv->NY+2)*3,MPI_FLOAT,gv->INDEX[1],gv->TAG1,gv->INDEX[2],gv->TAG1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&bufferrig_to_lef[0][1],(gv->NY+2)*3,MPI_FLOAT,gv->INDEX[2],gv->TAG2,gv->INDEX[1],gv->TAG2,MPI_COMM_WORLD,&status);

		

	/* no exchange if periodic boundary condition is applied */
	if ((gv->BOUNDARY) || (gv->POS[1]!=gv->NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=0;j<=gv->NY+1;j++){

			vx[j][gv->NX+1] = bufferlef_to_rig[j][1];
			vy[j][gv->NX+1] = bufferlef_to_rig[j][2];
			vz[j][gv->NX+1] = bufferlef_to_rig[j][3];
		
	}

						/* no exchange if periodic boundary condition is applied */
	if ((gv->BOUNDARY) || (gv->POS[1]!=0))	/* no boundary exchange at left edge of global grid */
	for (j=0;j<=gv->NY+1;j++){
			vx[j][0] = bufferrig_to_lef[j][1];
			vy[j][0] = bufferrig_to_lef[j][2];
			vz[j][0] = bufferrig_to_lef[j][3];
		
	}		
		
		
		
		
		
		
		
		
		
		
		
			
}
