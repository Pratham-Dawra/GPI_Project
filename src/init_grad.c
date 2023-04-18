
#include "fd.h"

void init_grad(int iter, MemInv * minv, GlobVar *gv, GlobVarInv *vinv)
{
    if (vinv->GRAD_METHOD == 2) {
        /* detect a change in inversion process and restart L-BFGS */
        if (iter == vinv->INV_RHO_ITER || iter == vinv->INV_VP_ITER || iter == vinv->INV_VS_ITER) {
            vinv->LBFGS_ITER_START = iter;
            if (vinv->WOLFE_CONDITION) {
                /* Restart Step Length search */
                vinv->ALPHA_SL_OLD = 1;
            }
            /* set values */
            vinv->FWI_RUN = 1;
            vinv->GRADIENT_OPTIMIZATION = 1;
        }

        /* restart L-BFGS */
        if (iter == vinv->LBFGS_ITER_START) {
            lbfgs_reset(iter, minv, gv, vinv);
            /* set values */
            vinv->FWI_RUN = 1;
            vinv->GRADIENT_OPTIMIZATION = 1;
        }
        /* Reset fail status of parabolic step length search */
        vinv->STEP3 = 0;
    }

    vinv->COUNTSTEP = 0;

    if (vinv->GRAD_METHOD == 1) {
        vinv->FWI_RUN = 1;
        vinv->STEPLENGTH_SEARCH = 0;
        vinv->GRADIENT_OPTIMIZATION = 1;
    }

}
