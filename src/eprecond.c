/*                          Acknowledgement                         */
/* This function is copied from DENISE Black Edition from D. Koehn  */
/* Licence: GNU GENERAL PUBLIC LICENSE Version 2, June 1991         */
/* https://github.com/daniel-koehn/DENISE-Black-Edition             */


#include "fd.h"

void eprecond(MemWavefield *mpw, float **W, GlobVar *gv)
{
    int i, j;

    for (j=1;j<=gv->NY;j++){
        for (i=1;i<=gv->NX;i++){
            W[j][i]+=((mpw->pvx[j][i] * mpw->pvx[j][i]) + (mpw->pvy[j][i] * mpw->pvy[j][i]));
        }
    }
}
