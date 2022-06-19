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
/* ------------------------------------------------------------------------
 * prepare of update of the stress tensor
 * ------------------------------------------------------------------------*/

#include "fd.h"

void prepare_update_s_vti(float *peta, float ** pc11, float **pc13, float ** pc33, float **pc55ipjp,
                          float **ptau11, float **ptau13, float ** ptau33, float **ptau55ipjp,
                          float ** pc55ipjpu, float ** pc13u, float **pc11u, float **pc33u,
                          float *** pc55ipjpd, float *** pc13d, float ***pc11d, float ***pc33d,
                          float *bip, float *cip) {

	extern int NX, NY, L;
	extern float DT;
	int i, j, l;



	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
            /* unrelaxed moduli */
			pc55ipjpu[j][i] = pc55ipjp[j][i]*DT*(1.0+L*ptau55ipjp[j][i]);
			pc13u[j][i] = pc13[j][i]*DT*(1.0+L*ptau13[j][i]);
            pc11u[j][i] = pc11[j][i]*DT*(1.0+L*ptau11[j][i]);
            pc33u[j][i] = pc33[j][i]*DT*(1.0+L*ptau33[j][i]);
			for (l=1;l<=L;l++){
				bip[l] = 1.0/(1.0+(peta[l]*0.5));
				cip[l] = 1.0-(peta[l]*0.5);
                /* module defects for each relaxation mechanism */
				pc55ipjpd[j][i][l] = pc55ipjp[j][i]*peta[l]*ptau55ipjp[j][i];
				pc13d[j][i][l] = pc13[j][i]*peta[l]*ptau13[j][i];
                pc33d[j][i][l] = pc33[j][i]*peta[l]*ptau33[j][i];
                pc11d[j][i][l] = pc11[j][i]*peta[l]*ptau11[j][i];
			}
		}
	}

}
