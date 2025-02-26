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
 *   initialisation of repeated comunications. This may reduce the
 *   network overhead. Communication is started later
 *   with MPI_START(request)
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, 
MPI_Request *req_send, MPI_Request *req_rec, GlobVar *gv){

	int fdo2;

	/* comunication initialisation for persistent communication */

	fdo2 = 2*(gv->FDORDER/2);

	/* buffer arrays are copied into local buffers using buffered send (bsend),
	  Actually send is activated by MPI_START(request) within time loop.
	  MPI_BSEND and MPI_RECV (see below) are then non-blocking.
	*/
	MPI_Bsend_init(&bufferlef_to_rig[1][1],gv->NY*fdo2,MPI_FLOAT,gv->INDEX[1],gv->TAG1,MPI_COMM_WORLD,&req_send[0]);
	MPI_Bsend_init(&bufferrig_to_lef[1][1],gv->NY*fdo2-1,MPI_FLOAT,gv->INDEX[2],gv->TAG2,MPI_COMM_WORLD,&req_send[1]);
	MPI_Bsend_init(&buffertop_to_bot[1][1],gv->NX*fdo2,MPI_FLOAT,gv->INDEX[3],gv->TAG3,MPI_COMM_WORLD,&req_send[2]);
	MPI_Bsend_init(&bufferbot_to_top[1][1],gv->NX*fdo2-1,MPI_FLOAT,gv->INDEX[4],gv->TAG4,MPI_COMM_WORLD,&req_send[3]);

	/* initialising of receive of buffer arrays. Same arrays for send and receive
	   are used. Thus, before starting receive, it is necessary to check if
	   Bsend has copied data into local buffers, i.e. has completed.
	*/
	MPI_Recv_init(&bufferlef_to_rig[1][1],gv->NY*fdo2,MPI_FLOAT,gv->INDEX[2],gv->TAG1,MPI_COMM_WORLD,&req_rec[0]);
	MPI_Recv_init(&bufferrig_to_lef[1][1],gv->NY*fdo2-1,MPI_FLOAT,gv->INDEX[1],gv->TAG2,MPI_COMM_WORLD,&req_rec[1]);
	MPI_Recv_init(&buffertop_to_bot[1][1],gv->NX*fdo2,MPI_FLOAT,gv->INDEX[4],gv->TAG3,MPI_COMM_WORLD,&req_rec[2]);
	MPI_Recv_init(&bufferbot_to_top[1][1],gv->NX*fdo2-1,MPI_FLOAT,gv->INDEX[3],gv->TAG4,MPI_COMM_WORLD,&req_rec[3]);

}
