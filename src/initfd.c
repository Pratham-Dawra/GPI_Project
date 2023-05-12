
#include "fd.h"
#include "logging.h"

static float *hc = NULL;

static float dhi = 0.f;

void op_s_fd2(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] = hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) * dhi;
    mpw->pvyx[j][i] = hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) * dhi;
    mpw->pvxy[j][i] = hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) * dhi;
    mpw->pvyy[j][i] = hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) * dhi;
}

void op_s_fd4(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] = (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2])) * dhi;
    mpw->pvyx[j][i] = (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1])) * dhi;
    mpw->pvxy[j][i] = (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i])) * dhi;
    mpw->pvyy[j][i] = (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i])) * dhi;
}

void op_s_fd6(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] =
        (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]) +
         hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3])) * dhi;
    mpw->pvyx[j][i] =
        (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]) +
         hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2])) * dhi;
    mpw->pvxy[j][i] =
        (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]) +
         hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i])) * dhi;
    mpw->pvyy[j][i] =
        (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]) +
         hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i])) * dhi;
}

void op_s_fd8(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] =
        (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]) +
         hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3]) + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4])) * dhi;
    mpw->pvyx[j][i] =
        (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]) +
         hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2]) + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3])) * dhi;
    mpw->pvxy[j][i] =
        (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]) +
         hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i]) + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i])) * dhi;
    mpw->pvyy[j][i] =
        (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]) +
         hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i]) + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i])) * dhi;
}

void op_s_fd10(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] =
        (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]) +
         hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3]) + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4]) +
         hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5])) * dhi;
    mpw->pvyy[j][i] =
        (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]) +
         hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i]) + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i]) +
         hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i])) * dhi;
    mpw->pvyx[j][i] =
        (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]) +
         hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2]) + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3]) +
         hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4])) * dhi;
    mpw->pvxy[j][i] =
        (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]) +
         hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i]) + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i]) +
         hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i])) * dhi;
}

void op_s_fd12(int i, int j, MemWavefield * mpw)
{
    mpw->pvxx[j][i] =
        (hc[1] * (mpw->pvx[j][i] - mpw->pvx[j][i - 1]) + hc[2] * (mpw->pvx[j][i + 1] - mpw->pvx[j][i - 2]) +
         hc[3] * (mpw->pvx[j][i + 2] - mpw->pvx[j][i - 3]) + hc[4] * (mpw->pvx[j][i + 3] - mpw->pvx[j][i - 4]) +
         hc[5] * (mpw->pvx[j][i + 4] - mpw->pvx[j][i - 5]) + hc[6] * (mpw->pvx[j][i + 5] - mpw->pvx[j][i - 6])) * dhi;
    mpw->pvyy[j][i] =
        (hc[1] * (mpw->pvy[j][i] - mpw->pvy[j - 1][i]) + hc[2] * (mpw->pvy[j + 1][i] - mpw->pvy[j - 2][i]) +
         hc[3] * (mpw->pvy[j + 2][i] - mpw->pvy[j - 3][i]) + hc[4] * (mpw->pvy[j + 3][i] - mpw->pvy[j - 4][i]) +
         hc[5] * (mpw->pvy[j + 4][i] - mpw->pvy[j - 5][i]) + hc[6] * (mpw->pvy[j + 5][i] - mpw->pvy[j - 6][i])) * dhi;
    mpw->pvyx[j][i] =
        (hc[1] * (mpw->pvy[j][i + 1] - mpw->pvy[j][i]) + hc[2] * (mpw->pvy[j][i + 2] - mpw->pvy[j][i - 1]) +
         hc[3] * (mpw->pvy[j][i + 3] - mpw->pvy[j][i - 2]) + hc[4] * (mpw->pvy[j][i + 4] - mpw->pvy[j][i - 3]) +
         hc[5] * (mpw->pvy[j][i + 5] - mpw->pvy[j][i - 4]) + hc[6] * (mpw->pvy[j][i + 6] - mpw->pvy[j][i - 5])) * dhi;
    mpw->pvxy[j][i] =
        (hc[1] * (mpw->pvx[j + 1][i] - mpw->pvx[j][i]) + hc[2] * (mpw->pvx[j + 2][i] - mpw->pvx[j - 1][i]) +
         hc[3] * (mpw->pvx[j + 3][i] - mpw->pvx[j - 2][i]) + hc[4] * (mpw->pvx[j + 4][i] - mpw->pvx[j - 3][i]) +
         hc[5] * (mpw->pvx[j + 5][i] - mpw->pvx[j - 4][i]) + hc[6] * (mpw->pvx[j + 6][i] - mpw->pvx[j - 5][i])) * dhi;
}

void op_v_fd2(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]);
    *sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]);
    *sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]);
    *syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]);
}

void op_v_fd4(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x = hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]) + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]);
    *sxy_x = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]) + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]);
    *sxy_y = hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]) + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]);
    *syy_y = hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]) + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]);
}

void op_v_fd6(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x =
        hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]) + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]) +
        hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2]);
    *sxy_x =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]) + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]) +
        hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3]);
    *sxy_y =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]) + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]) +
        hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i]);
    *syy_y =
        hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]) + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]) +
        hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i]);
}

void op_v_fd8(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x =
        hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]) + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]) +
        hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2]) + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3]);
    *sxy_x =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]) + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]) +
        hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3]) + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4]);
    *sxy_y =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]) + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]) +
        hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i]) + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i]);
    *syy_y =
        hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]) + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]) +
        hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i]) + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i]);
}

void op_v_fd10(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x =
        hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]) + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]) +
        hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2]) + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3]) +
        hc[5] * (mpw->psxx[j][i + 5] - mpw->psxx[j][i - 4]);
    *sxy_x =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]) + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]) +
        hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3]) + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4]) +
        hc[5] * (mpw->psxy[j][i + 4] - mpw->psxy[j][i - 5]);
    *sxy_y =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]) + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]) +
        hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i]) + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i]) +
        hc[5] * (mpw->psxy[j + 4][i] - mpw->psxy[j - 5][i]);
    *syy_y =
        hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]) + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]) +
        hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i]) + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i]) +
        hc[5] * (mpw->psyy[j + 5][i] - mpw->psyy[j - 4][i]);
}

void op_v_fd12(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield * mpw)
{
    *sxx_x =
        hc[1] * (mpw->psxx[j][i + 1] - mpw->psxx[j][i]) + hc[2] * (mpw->psxx[j][i + 2] - mpw->psxx[j][i - 1]) +
        hc[3] * (mpw->psxx[j][i + 3] - mpw->psxx[j][i - 2]) + hc[4] * (mpw->psxx[j][i + 4] - mpw->psxx[j][i - 3]) +
        hc[5] * (mpw->psxx[j][i + 5] - mpw->psxx[j][i - 4]) + hc[6] * (mpw->psxx[j][i + 6] - mpw->psxx[j][i - 5]);
    *sxy_x =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j][i - 1]) + hc[2] * (mpw->psxy[j][i + 1] - mpw->psxy[j][i - 2]) +
        hc[3] * (mpw->psxy[j][i + 2] - mpw->psxy[j][i - 3]) + hc[4] * (mpw->psxy[j][i + 3] - mpw->psxy[j][i - 4]) +
        hc[5] * (mpw->psxy[j][i + 4] - mpw->psxy[j][i - 5]) + hc[6] * (mpw->psxy[j][i + 5] - mpw->psxy[j][i - 6]);
    *sxy_y =
        hc[1] * (mpw->psxy[j][i] - mpw->psxy[j - 1][i]) + hc[2] * (mpw->psxy[j + 1][i] - mpw->psxy[j - 2][i]) +
        hc[3] * (mpw->psxy[j + 2][i] - mpw->psxy[j - 3][i]) + hc[4] * (mpw->psxy[j + 3][i] - mpw->psxy[j - 4][i]) +
        hc[5] * (mpw->psxy[j + 4][i] - mpw->psxy[j - 5][i]) + hc[6] * (mpw->psxy[j + 5][i] - mpw->psxy[j - 6][i]);
    *syy_y =
        hc[1] * (mpw->psyy[j + 1][i] - mpw->psyy[j][i]) + hc[2] * (mpw->psyy[j + 2][i] - mpw->psyy[j - 1][i]) +
        hc[3] * (mpw->psyy[j + 3][i] - mpw->psyy[j - 2][i]) + hc[4] * (mpw->psyy[j + 4][i] - mpw->psyy[j - 3][i]) +
        hc[5] * (mpw->psyy[j + 5][i] - mpw->psyy[j - 4][i]) + hc[6] * (mpw->psyy[j + 6][i] - mpw->psyy[j - 5][i]);
}

static void update_fd_fct_ptr(int order, GlobVar * gv)
{
    switch (order) {
      case 2:
          gv->FDOP_S = &op_s_fd2;
          gv->FDOP_V = &op_v_fd2;
          break;
      case 4:
          gv->FDOP_S = &op_s_fd4;
          gv->FDOP_V = &op_v_fd4;
          break;
      case 6:
          gv->FDOP_S = &op_s_fd6;
          gv->FDOP_V = &op_v_fd6;
          break;
      case 8:
          gv->FDOP_S = &op_s_fd8;
          gv->FDOP_V = &op_v_fd8;
          break;
      case 10:
          gv->FDOP_S = &op_s_fd10;
          gv->FDOP_V = &op_v_fd10;
          break;
      case 12:
          gv->FDOP_S = &op_s_fd12;
          gv->FDOP_V = &op_v_fd12;
          break;
      default:
          log_fatal("Unsupported FDORDER encountered.\n");
          break;
    }
    return;
}

void initfd(GlobVar * gv)
{
    hc = holbergcoeff(gv);
    dhi = 1.f / gv->DH;

    update_fd_fct_ptr(gv->FDORDER, gv);

    return;
}

void set_fd_order(int new_order, GlobVar * gv)
{
    gv->FDORDER = new_order;

    update_fd_fct_ptr(new_order, gv);

    return;
}

int get_fd_order(GlobVar * gv)
{
    return gv->FDORDER;
}
