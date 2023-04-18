
#include "fd.h"

void init_inv(int iter, MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{
    int i, j;

    /* initialization of L2 calculation */
    vinv->L2 = 0.0;
    vinv->ENERGY = 0.0;
    vinv->L2_ALL_SHOTS = 0.0;
    vinv->ENERGY_ALL_SHOTS = 0.0;
    vinv->KILLED_TRACES = 0;
    vinv->KILLED_TRACES_TESTSHOTS = 0;

    vinv->EPSILON = 0.0;        /* test step length */
    exchange_par(gv, vinv);     // ??? necessary

    /* initialize waveconv matrix */
    for (j = 1; j <= gv->NY; j++) {
        for (i = 1; i <= gv->NX; i++) {
            minv->waveconv[j][i] = 0.0;
            minv->waveconv_rho[j][i] = 0.0;
            minv->waveconv_u[j][i] = 0.0;
        }
    }

    if ((vinv->EPRECOND > 0) && (vinv->EPRECOND_ITER == iter || (vinv->EPRECOND_ITER == 0))) {
        for (j = 1; j <= gv->NY; j++) {
            for (i = 1; i <= gv->NX; i++) {
                minv->We_sum[j][i] = 0.0;
            }
        }
    }

    vinv->ITESTSHOT = vinv->TESTSHOT_START;
    vinv->SWS_TESTSHOT = 0;
}
