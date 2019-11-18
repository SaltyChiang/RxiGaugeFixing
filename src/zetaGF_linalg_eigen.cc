#include <iostream>

#include "../include/zetaGF_lattice.h"
#include "../include/zetaGF_linalg.h"
#include "../extension/Eigen/Dense"

void zgfGenAField_eigen(zgfGaugeMatrix *af, zgfGaugeMatrix *gf)
{
  Eigen::Matrix3cd *afEigen = (Eigen::Matrix3cd *)af;
  Eigen::Matrix3cd *gfEigen = (Eigen::Matrix3cd *)gf;
  Eigen::dcomplex coef(0.0, 2.0);
#pragma omp parallel for
  for (int i = 0; i < VOL * Nd; i++)
  {
    afEigen[i] = (gfEigen[i] - gfEigen[i].adjoint()) / coef;
    afEigen[i] -= afEigen[i].trace() / 3.0 * Eigen::Matrix3cd::Identity();
  }
}

void zgfGenDeltaField_eigen(zgfGaugeMatrix *df, zgfGaugeMatrix *af)
{
  Eigen::Matrix3cd *dfEigen = (Eigen::Matrix3cd *)df;
  Eigen::Matrix3cd *afEigen = (Eigen::Matrix3cd *)af;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    int t = i / (Nz * Ny * Nz);
    int z = i % (Nz * Ny * Nx) / (Ny * Nx);
    int y = i % (Ny * Nx) / (Nx);
    int x = i % (Nx);
    int tm = i + (t == 0 ? (Nt - 1) : -1) * Nz * Ny * Nx;
    int zm = i + (z == 0 ? (Nz - 1) : -1) * Ny * Nx;
    int ym = i + (y == 0 ? (Ny - 1) : -1) * Nx;
    int xm = i + (x == 0 ? (Nx - 1) : -1);

    dfEigen[i] = afEigen[i * Nd + 0] - afEigen[xm * Nd + 0] +
                 afEigen[i * Nd + 1] - afEigen[ym * Nd + 1] +
                 afEigen[i * Nd + 2] - afEigen[zm * Nd + 2] +
                 afEigen[i * Nd + 3] - afEigen[tm * Nd + 3];
  }
}

void zgfInitGaugeRotateField_eigen(zgfGaugeMatrix *grf)
{
  Eigen::Matrix3cd *grfEigen = (Eigen::Matrix3cd *)grf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
    grfEigen[i] = Eigen::Matrix3cd::Identity();
}

void zgfGenKField_eigen(zgfGaugeMatrix *kf, zgfGaugeMatrix *gf, zgfGaugeMatrix *grf)
{
  Eigen::Matrix3cd *kfEigen = (Eigen::Matrix3cd *)kf;
  Eigen::Matrix3cd *gfEigen = (Eigen::Matrix3cd *)gf;
  Eigen::Matrix3cd *grfEigen = (Eigen::Matrix3cd *)grf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    int t = i / (Nz * Ny * Nz);
    int z = i % (Nz * Ny * Nx) / (Ny * Nx);
    int y = i % (Ny * Nx) / (Nx);
    int x = i % (Nx);
    int tm = i + (t == 0 ? (Nt - 1) : -1) * Nz * Ny * Nx;
    int zm = i + (z == 0 ? (Nz - 1) : -1) * Ny * Nx;
    int ym = i + (y == 0 ? (Ny - 1) : -1) * Nx;
    int xm = i + (x == 0 ? (Nx - 1) : -1);
    int tp = i + (t == (Nt - 1) ? (1 - Nt) : 1) * Nz * Ny * Nx;
    int zp = i + (z == (Nz - 1) ? (1 - Nz) : 1) * Ny * Nx;
    int yp = i + (y == (Ny - 1) ? (1 - Ny) : 1) * Nx;
    int xp = i + (x == (Nx - 1) ? (1 - Nx) : 1);

    kfEigen[i] = grfEigen[xp] * gfEigen[i * Nd + 0].adjoint() +
                 grfEigen[xm] * gfEigen[xm * Nd + 0] +
                 grfEigen[yp] * gfEigen[i * Nd + 1].adjoint() +
                 grfEigen[ym] * gfEigen[ym * Nd + 1] +
                 grfEigen[zp] * gfEigen[i * Nd + 2].adjoint() +
                 grfEigen[zm] * gfEigen[zm * Nd + 2] +
                 grfEigen[tp] * gfEigen[i * Nd + 3].adjoint() +
                 grfEigen[tm] * gfEigen[tm * Nd + 3];
  }
}

void zgfGenGaugeRotateField_eigen(zgfGaugeMatrix *grf, zgfGaugeMatrix *kf)
{
  Eigen::Matrix3cd *grfEigen = (Eigen::Matrix3cd *)grf;
  Eigen::Matrix3cd *kfEigen = (Eigen::Matrix3cd *)kf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    grfEigen[i] = kfEigen[i] / sqrt(kfEigen[i].determinant());
    grfEigen[i] = grfEigen[i] * grfEigen[i];
  }
}

void zgfUpdateGaugeField_eigen(zgfGaugeMatrix *gf, zgfGaugeMatrix *grf)
{
  Eigen::Matrix3cd *gfEigen = (Eigen::Matrix3cd *)gf;
  Eigen::Matrix3cd *grfEigen = (Eigen::Matrix3cd *)grf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    int t = i / (Nz * Ny * Nz);
    int z = i % (Nz * Ny * Nx) / (Ny * Nx);
    int y = i % (Ny * Nx) / (Nx);
    int x = i % (Nx);
    // int tm = i + (t == 0 ? (Nt - 1) : -1) * Nz * Ny * Nx;
    // int zm = i + (z == 0 ? (Nz - 1) : -1) * Ny * Nx;
    // int ym = i + (y == 0 ? (Ny - 1) : -1) * Nx;
    // int xm = i + (x == 0 ? (Nx - 1) : -1);
    int tp = i + (t == (Nt - 1) ? (1 - Nt) : 1) * Nz * Ny * Nx;
    int zp = i + (z == (Nz - 1) ? (1 - Nz) : 1) * Ny * Nx;
    int yp = i + (y == (Ny - 1) ? (1 - Ny) : 1) * Nx;
    int xp = i + (x == (Nx - 1) ? (1 - Nx) : 1);

    gfEigen[i * Nd + 0] = grfEigen[i] * gfEigen[i * Nd + 0] * grfEigen[xp].adjoint();
    gfEigen[i * Nd + 1] = grfEigen[i] * gfEigen[i * Nd + 1] * grfEigen[yp].adjoint();
    gfEigen[i * Nd + 2] = grfEigen[i] * gfEigen[i * Nd + 2] * grfEigen[zp].adjoint();
    gfEigen[i * Nd + 3] = grfEigen[i] * gfEigen[i * Nd + 3] * grfEigen[tp].adjoint();
  }
}

double zgfGetTheta_eigen(zgfGaugeMatrix *df)
{
  Eigen::Matrix3cd *dfEigen = (Eigen::Matrix3cd *)df;
  double *theta;
  theta = (double *)malloc(VOL * sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
    theta[i] = dfEigen[i].squaredNorm();
  for (int i = 1; i < VOL; i++)
    theta[0] += theta[i];
  return theta[0] / Nc / VOL;
}
