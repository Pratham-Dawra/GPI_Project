
#include "fd.h"
#include "logging.h"

static float *hc = NULL;

static float dhi = 0.f;

void op_s_fd2(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx = hc[1] * (vx[j][i] - vx[j][i - 1]) * dhi;
    *vyx = hc[1] * (vy[j][i + 1] - vy[j][i]) * dhi;
    *vxy = hc[1] * (vx[j + 1][i] - vx[j][i]) * dhi;
    *vyy = hc[1] * (vy[j][i] - vy[j - 1][i]) * dhi;
}

void op_s_fd4(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx = (hc[1] * (vx[j][i] - vx[j][i - 1]) + hc[2] * (vx[j][i + 1] - vx[j][i - 2])) * dhi;
    *vyx = (hc[1] * (vy[j][i + 1] - vy[j][i]) + hc[2] * (vy[j][i + 2] - vy[j][i - 1])) * dhi;
    *vxy = (hc[1] * (vx[j + 1][i] - vx[j][i]) + hc[2] * (vx[j + 2][i] - vx[j - 1][i])) * dhi;
    *vyy = (hc[1] * (vy[j][i] - vy[j - 1][i]) + hc[2] * (vy[j + 1][i] - vy[j - 2][i])) * dhi;
}

void op_s_fd6(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx =
        (hc[1] * (vx[j][i] - vx[j][i - 1]) + hc[2] * (vx[j][i + 1] - vx[j][i - 2]) +
         hc[3] * (vx[j][i + 2] - vx[j][i - 3])) * dhi;
    *vyx =
        (hc[1] * (vy[j][i + 1] - vy[j][i]) + hc[2] * (vy[j][i + 2] - vy[j][i - 1]) +
         hc[3] * (vy[j][i + 3] - vy[j][i - 2])) * dhi;
    *vxy =
        (hc[1] * (vx[j + 1][i] - vx[j][i]) + hc[2] * (vx[j + 2][i] - vx[j - 1][i]) +
         hc[3] * (vx[j + 3][i] - vx[j - 2][i])) * dhi;
    *vyy =
        (hc[1] * (vy[j][i] - vy[j - 1][i]) + hc[2] * (vy[j + 1][i] - vy[j - 2][i]) +
         hc[3] * (vy[j + 2][i] - vy[j - 3][i])) * dhi;
}

void op_s_fd8(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx =
        (hc[1] * (vx[j][i] - vx[j][i - 1]) + hc[2] * (vx[j][i + 1] - vx[j][i - 2]) +
         hc[3] * (vx[j][i + 2] - vx[j][i - 3]) + hc[4] * (vx[j][i + 3] - vx[j][i - 4])) * dhi;
    *vyx =
        (hc[1] * (vy[j][i + 1] - vy[j][i]) + hc[2] * (vy[j][i + 2] - vy[j][i - 1]) +
         hc[3] * (vy[j][i + 3] - vy[j][i - 2]) + hc[4] * (vy[j][i + 4] - vy[j][i - 3])) * dhi;
    *vxy =
        (hc[1] * (vx[j + 1][i] - vx[j][i]) + hc[2] * (vx[j + 2][i] - vx[j - 1][i]) +
         hc[3] * (vx[j + 3][i] - vx[j - 2][i]) + hc[4] * (vx[j + 4][i] - vx[j - 3][i])) * dhi;
    *vyy =
        (hc[1] * (vy[j][i] - vy[j - 1][i]) + hc[2] * (vy[j + 1][i] - vy[j - 2][i]) +
         hc[3] * (vy[j + 2][i] - vy[j - 3][i]) + hc[4] * (vy[j + 3][i] - vy[j - 4][i])) * dhi;
}

void op_s_fd10(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx =
        (hc[1] * (vx[j][i] - vx[j][i - 1]) + hc[2] * (vx[j][i + 1] - vx[j][i - 2]) +
         hc[3] * (vx[j][i + 2] - vx[j][i - 3]) + hc[4] * (vx[j][i + 3] - vx[j][i - 4]) + hc[5] * (vx[j][i + 4] -
                                                                                                  vx[j][i - 5])) * dhi;
    *vyy =
        (hc[1] * (vy[j][i] - vy[j - 1][i]) + hc[2] * (vy[j + 1][i] - vy[j - 2][i]) +
         hc[3] * (vy[j + 2][i] - vy[j - 3][i]) + hc[4] * (vy[j + 3][i] - vy[j - 4][i]) + hc[5] * (vy[j + 4][i] -
                                                                                                  vy[j - 5][i])) * dhi;
    *vyx =
        (hc[1] * (vy[j][i + 1] - vy[j][i]) + hc[2] * (vy[j][i + 2] - vy[j][i - 1]) +
         hc[3] * (vy[j][i + 3] - vy[j][i - 2]) + hc[4] * (vy[j][i + 4] - vy[j][i - 3]) + hc[5] * (vy[j][i + 5] -
                                                                                                  vy[j][i - 4])) * dhi;
    *vxy =
        (hc[1] * (vx[j + 1][i] - vx[j][i]) + hc[2] * (vx[j + 2][i] - vx[j - 1][i]) +
         hc[3] * (vx[j + 3][i] - vx[j - 2][i]) + hc[4] * (vx[j + 4][i] - vx[j - 3][i]) + hc[5] * (vx[j + 5][i] -
                                                                                                  vx[j - 4][i])) * dhi;
}

void op_s_fd12(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, float **vx, float **vy)
{
    *vxx =
        (hc[1] * (vx[j][i] - vx[j][i - 1]) + hc[2] * (vx[j][i + 1] - vx[j][i - 2]) +
         hc[3] * (vx[j][i + 2] - vx[j][i - 3]) + hc[4] * (vx[j][i + 3] - vx[j][i - 4]) + hc[5] * (vx[j][i + 4] -
                                                                                                  vx[j][i - 5]) +
         hc[6] * (vx[j][i + 5] - vx[j][i - 6])) * dhi;
    *vyy =
        (hc[1] * (vy[j][i] - vy[j - 1][i]) + hc[2] * (vy[j + 1][i] - vy[j - 2][i]) +
         hc[3] * (vy[j + 2][i] - vy[j - 3][i]) + hc[4] * (vy[j + 3][i] - vy[j - 4][i]) + hc[5] * (vy[j + 4][i] -
                                                                                                  vy[j - 5][i]) +
         hc[6] * (vy[j + 5][i] - vy[j - 6][i])) * dhi;
    *vyx =
        (hc[1] * (vy[j][i + 1] - vy[j][i]) + hc[2] * (vy[j][i + 2] - vy[j][i - 1]) +
         hc[3] * (vy[j][i + 3] - vy[j][i - 2]) + hc[4] * (vy[j][i + 4] - vy[j][i - 3]) + hc[5] * (vy[j][i + 5] -
                                                                                                  vy[j][i - 4]) +
         hc[6] * (vy[j][i + 6] - vy[j][i - 5])) * dhi;
    *vxy =
        (hc[1] * (vx[j + 1][i] - vx[j][i]) + hc[2] * (vx[j + 2][i] - vx[j - 1][i]) +
         hc[3] * (vx[j + 3][i] - vx[j - 2][i]) + hc[4] * (vx[j + 4][i] - vx[j - 3][i]) + hc[5] * (vx[j + 5][i] -
                                                                                                  vx[j - 4][i]) +
         hc[6] * (vx[j + 6][i] - vx[j - 5][i])) * dhi;
}

void op_v_fd2(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
              float **sxy)
{
    *sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i]);
    *sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1]);
    *sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i]);
    *syy_y = hc[1] * (syy[j + 1][i] - syy[j][i]);
}

void op_v_fd4(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
              float **sxy)
{
    *sxx_x = hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]);
    *sxy_x = hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]);
    *sxy_y = hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]);
    *syy_y = hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]);
}

void op_v_fd6(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
              float **sxy)
{
    *sxx_x =
        hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]) + hc[3] * (sxx[j][i + 3] -
                                                                                                 sxx[j][i - 2]);
    *sxy_x =
        hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]) + hc[3] * (sxy[j][i + 2] -
                                                                                                 sxy[j][i - 3]);
    *sxy_y =
        hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]) + hc[3] * (sxy[j + 2][i] -
                                                                                                 sxy[j - 3][i]);
    *syy_y =
        hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]) + hc[3] * (syy[j + 3][i] -
                                                                                                 syy[j - 2][i]);
}

void op_v_fd8(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
              float **sxy)
{
    *sxx_x =
        hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]) + hc[3] * (sxx[j][i + 3] -
                                                                                                 sxx[j][i - 2]) +
        hc[4] * (sxx[j][i + 4] - sxx[j][i - 3]);
    *sxy_x =
        hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]) + hc[3] * (sxy[j][i + 2] -
                                                                                                 sxy[j][i - 3]) +
        hc[4] * (sxy[j][i + 3] - sxy[j][i - 4]);
    *sxy_y =
        hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]) + hc[3] * (sxy[j + 2][i] -
                                                                                                 sxy[j - 3][i]) +
        hc[4] * (sxy[j + 3][i] - sxy[j - 4][i]);
    *syy_y =
        hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]) + hc[3] * (syy[j + 3][i] -
                                                                                                 syy[j - 2][i]) +
        hc[4] * (syy[j + 4][i] - syy[j - 3][i]);
}

void op_v_fd10(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
               float **sxy)
{
    *sxx_x =
        hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]) + hc[3] * (sxx[j][i + 3] -
                                                                                                 sxx[j][i - 2]) +
        hc[4] * (sxx[j][i + 4] - sxx[j][i - 3]) + hc[5] * (sxx[j][i + 5] - sxx[j][i - 4]);
    *sxy_x =
        hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]) + hc[3] * (sxy[j][i + 2] -
                                                                                                 sxy[j][i - 3]) +
        hc[4] * (sxy[j][i + 3] - sxy[j][i - 4]) + hc[5] * (sxy[j][i + 4] - sxy[j][i - 5]);
    *sxy_y =
        hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]) + hc[3] * (sxy[j + 2][i] -
                                                                                                 sxy[j - 3][i]) +
        hc[4] * (sxy[j + 3][i] - sxy[j - 4][i]) + hc[5] * (sxy[j + 4][i] - sxy[j - 5][i]);
    *syy_y =
        hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]) + hc[3] * (syy[j + 3][i] -
                                                                                                 syy[j - 2][i]) +
        hc[4] * (syy[j + 4][i] - syy[j - 3][i]) + hc[5] * (syy[j + 5][i] - syy[j - 4][i]);
}

void op_v_fd12(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, float **sxx, float **syy,
               float **sxy)
{
    *sxx_x =
        hc[1] * (sxx[j][i + 1] - sxx[j][i]) + hc[2] * (sxx[j][i + 2] - sxx[j][i - 1]) + hc[3] * (sxx[j][i + 3] -
                                                                                                 sxx[j][i - 2]) +
        hc[4] * (sxx[j][i + 4] - sxx[j][i - 3]) + hc[5] * (sxx[j][i + 5] - sxx[j][i - 4]) + hc[6] * (sxx[j][i + 6] -
                                                                                                     sxx[j][i - 5]);
    *sxy_x =
        hc[1] * (sxy[j][i] - sxy[j][i - 1]) + hc[2] * (sxy[j][i + 1] - sxy[j][i - 2]) + hc[3] * (sxy[j][i + 2] -
                                                                                                 sxy[j][i - 3]) +
        hc[4] * (sxy[j][i + 3] - sxy[j][i - 4]) + hc[5] * (sxy[j][i + 4] - sxy[j][i - 5]) + hc[6] * (sxy[j][i + 5] -
                                                                                                     sxy[j][i - 6]);
    *sxy_y =
        hc[1] * (sxy[j][i] - sxy[j - 1][i]) + hc[2] * (sxy[j + 1][i] - sxy[j - 2][i]) + hc[3] * (sxy[j + 2][i] -
                                                                                                 sxy[j - 3][i]) +
        hc[4] * (sxy[j + 3][i] - sxy[j - 4][i]) + hc[5] * (sxy[j + 4][i] - sxy[j - 5][i]) + hc[6] * (sxy[j + 5][i] -
                                                                                                     sxy[j - 6][i]);
    *syy_y =
        hc[1] * (syy[j + 1][i] - syy[j][i]) + hc[2] * (syy[j + 2][i] - syy[j - 1][i]) + hc[3] * (syy[j + 3][i] -
                                                                                                 syy[j - 2][i]) +
        hc[4] * (syy[j + 4][i] - syy[j - 3][i]) + hc[5] * (syy[j + 5][i] - syy[j - 4][i]) + hc[6] * (syy[j + 6][i] -
                                                                                                     syy[j - 5][i]);
}

static void update_fd_fct_ptr(int order, GlobVar *gv)
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

void initfd(GlobVar *gv)
{
    hc = holbergcoeff(gv);
    dhi = 1.f / gv->DH;

    update_fd_fct_ptr(gv->FDORDER, gv);

    return;
}

void set_fd_order(int new_order, GlobVar *gv)
{
    gv->FDORDER = new_order;

    update_fd_fct_ptr(new_order, gv);

    return;
}

int get_fd_order(GlobVar *gv)
{
    return gv->FDORDER;
}
