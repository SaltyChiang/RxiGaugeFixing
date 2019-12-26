#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_linalg.h"
#include "../extension/Eigen/Dense"

namespace RxiGF
{

typedef Eigen::Matrix<Eigen::dcomplex, Nc, Nc, Eigen::RowMajor> ColorMatrix_Eigen;
typedef Eigen::Matrix<Eigen::dcomplex, Nc, 1> ColorVector_Eigen;

void CopyGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf)
{
  memcpy((void *)tgf, (void *)gf, sizeof(ColorMatrix) * VOL * Nd);
}

void InitGaugeRotateField_eigen(ColorMatrix *grf)
{
  ColorMatrix_Eigen *grfEigen = (ColorMatrix_Eigen *)grf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
    grfEigen[i] = ColorMatrix_Eigen::Identity();
}

void ReunitGaugeRotateField_eigen(ColorMatrix *grf)
{
  ColorMatrix_Eigen *grfEigen = (ColorMatrix_Eigen *)grf;
  ColorVector_Eigen *u = (ColorVector_Eigen *)malloc(sizeof(ColorVector_Eigen) * VOL);
  ColorVector_Eigen *v = (ColorVector_Eigen *)malloc(sizeof(ColorVector_Eigen) * VOL);
  ColorVector_Eigen *w = (ColorVector_Eigen *)malloc(sizeof(ColorVector_Eigen) * VOL);
  double *t1 = (double *)malloc(sizeof(double) * VOL);
  Eigen::dcomplex *t2 = (Eigen::dcomplex *)malloc(sizeof(Eigen::dcomplex) * VOL);
  double *t3 = (double *)malloc(sizeof(double) * VOL);
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    u[i] << grfEigen[i](0, 0), grfEigen[i](1, 0), grfEigen[i](2, 0);
    v[i] << grfEigen[i](0, 1), grfEigen[i](1, 1), grfEigen[i](2, 1);
    t1[i] = u[i].norm();
    u[i] /= t1[i];
    t2[i] = u[i].dot(v[i]);
    v[i] -= t2[i] * u[i];
    t3[i] = v[i].norm();
    v[i] /= t3[i];
    grfEigen[i](0, 0) = u[i](0);
    grfEigen[i](1, 0) = u[i](1);
    grfEigen[i](2, 0) = u[i](2);
    grfEigen[i](0, 1) = v[i](0);
    grfEigen[i](1, 1) = v[i](1);
    grfEigen[i](2, 1) = v[i](2);
    w[i](0) = u[i](1) * v[i](2) - u[i](2) * v[i](1);
    w[i](1) = u[i](2) * v[i](0) - u[i](0) * v[i](2);
    w[i](2) = u[i](0) * v[i](1) - u[i](1) * v[i](0);
    w[i] = w[i].conjugate();
    grfEigen[i](0, 2) = w[i](0);
    grfEigen[i](1, 2) = w[i](1);
    grfEigen[i](2, 2) = w[i](2);
  }
  free(u);
  free(v);
  free(w);
  free(t1);
  free(t2);
  free(t3);
}

void UpdateGaugeField_eigen(ColorMatrix *gf, ColorMatrix *grf)
{
  ColorMatrix_Eigen *gfEigen = (ColorMatrix_Eigen *)gf;
  ColorMatrix_Eigen *grfEigen = (ColorMatrix_Eigen *)grf;
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

void UpdateGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf, ColorMatrix *grf)
{

  ColorMatrix_Eigen *tgfEigen = (ColorMatrix_Eigen *)tgf;
  ColorMatrix_Eigen *gfEigen = (ColorMatrix_Eigen *)gf;
  ColorMatrix_Eigen *grfEigen = (ColorMatrix_Eigen *)grf;
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    int t = i / (Nz * Ny * Nz);
    int z = i % (Nz * Ny * Nx) / (Ny * Nx);
    int y = i % (Ny * Nx) / (Nx);
    int x = i % (Nx);
    int tp = i + (t == (Nt - 1) ? (1 - Nt) : 1) * Nz * Ny * Nx;
    int zp = i + (z == (Nz - 1) ? (1 - Nz) : 1) * Ny * Nx;
    int yp = i + (y == (Ny - 1) ? (1 - Ny) : 1) * Nx;
    int xp = i + (x == (Nx - 1) ? (1 - Nx) : 1);

    tgfEigen[i * Nd + 0] = grfEigen[i] * gfEigen[i * Nd + 0] * grfEigen[xp].adjoint();
    tgfEigen[i * Nd + 1] = grfEigen[i] * gfEigen[i * Nd + 1] * grfEigen[yp].adjoint();
    tgfEigen[i * Nd + 2] = grfEigen[i] * gfEigen[i * Nd + 2] * grfEigen[zp].adjoint();
    tgfEigen[i * Nd + 3] = grfEigen[i] * gfEigen[i * Nd + 3] * grfEigen[tp].adjoint();
  }
}

void GenAField_eigen(ColorMatrix *af, ColorMatrix *gf)
{
  ColorMatrix_Eigen *afEigen = (ColorMatrix_Eigen *)af;
  ColorMatrix_Eigen *gfEigen = (ColorMatrix_Eigen *)gf;
  Eigen::dcomplex coef(0.0, 2.0);
#pragma omp parallel for
  for (int i = 0; i < VOL * Nd; i++)
  {
    afEigen[i] = (gfEigen[i] - gfEigen[i].adjoint()) / coef;
    afEigen[i] -= afEigen[i].trace() / 3.0 * ColorMatrix_Eigen::Identity();
  }
}

void GenDeltaField_eigen(ColorMatrix *df, ColorMatrix *af)
{
  ColorMatrix_Eigen *dfEigen = (ColorMatrix_Eigen *)df;
  ColorMatrix_Eigen *afEigen = (ColorMatrix_Eigen *)af;
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

double GetTheta_eigen(ColorMatrix *df, ColorMatrix *af, ColorMatrix *gf)
{
  GenAField_eigen(af, gf);
  GenDeltaField_eigen(df, af);
  ColorMatrix_Eigen *dfEigen = (ColorMatrix_Eigen *)df;
  double *theta;
  theta = (double *)malloc(VOL * sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
    theta[i] = dfEigen[i].squaredNorm();
  for (int i = 1; i < VOL; i++)
    theta[0] += theta[i];
  double thetaOut = theta[0] / (VOL * Nc);
  free(theta);
  return thetaOut;
}

double GetFunctional_eigen(ColorMatrix *gf)
{
  ColorMatrix_Eigen *gfEigen = (ColorMatrix_Eigen *)gf;
  double *func;
  func = (double *)malloc(VOL * Nd * sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < VOL * Nd; i++)
    func[i] = gfEigen[i].trace().real();
  for (int i = 1; i < VOL * Nd; i++)
    func[0] += func[i];
  double funcOut = func[0] / (VOL * Nd * Nc);
  free(func);
  return funcOut;
}

/*
 * Deprecated.
 */
void GenKField_eigen(ColorMatrix *kf, ColorMatrix *gf, ColorMatrix *grf)
{
  ColorMatrix_Eigen *kfEigen = (ColorMatrix_Eigen *)kf;
  ColorMatrix_Eigen *gfEigen = (ColorMatrix_Eigen *)gf;
  ColorMatrix_Eigen *grfEigen = (ColorMatrix_Eigen *)grf;
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

} // namespace RxiGF