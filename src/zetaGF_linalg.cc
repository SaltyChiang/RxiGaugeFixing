#include <cstdlib>

#include "../include/zetaGF_lattice.h"
#include "../include/zetaGF_linalg.h"

zgfComplex operator+(zgfComplex a, zgfComplex b)
{
  zgfComplex c;
  c.re = a.re + b.re;
  c.im = a.im + b.im;
  return c;
}

zgfComplex operator-(zgfComplex a, zgfComplex b)
{
  zgfComplex c;
  c.re = a.re - b.re;
  c.im = a.im - b.im;
  return c;
}

zgfComplex operator*(zgfComplex a, zgfComplex b)
{
  zgfComplex c;
  c.re = a.re * b.re - a.im * b.im;
  c.im = a.re + b.im + a.im * b.re;
  return c;
}

zgfGaugeMatrix operator+(zgfGaugeMatrix a, zgfGaugeMatrix b)
{
  PREC *cFloat = (PREC *)malloc(sizeof(zgfGaugeMatrix));
  PREC *aFloat = (PREC *)&a;
  PREC *bFloat = (PREC *)&b;
  for (int i = 0; i < Nc * Nc * 2; i++)
    cFloat[i] = aFloat[i] + bFloat[i];
  return *(zgfGaugeMatrix *)cFloat;
}

zgfGaugeMatrix operator-(zgfGaugeMatrix a, zgfGaugeMatrix b)
{
  PREC *cFloat = (PREC *)malloc(sizeof(zgfGaugeMatrix));
  PREC *aFloat = (PREC *)&a;
  PREC *bFloat = (PREC *)&b;
  for (int i = 0; i < Nc * Nc * 2; i++)
    cFloat[i] = aFloat[i] - bFloat[i];
  return *(zgfGaugeMatrix *)cFloat;
}

zgfGaugeMatrix operator*(zgfGaugeMatrix a, zgfGaugeMatrix b)
{
  PREC *cFloat = (PREC *)malloc(sizeof(zgfGaugeMatrix));
  PREC *aFloat = (PREC *)&a;
  PREC *bFloat = (PREC *)&b;
  for (int i = 0; i < Nc * Nc * 2; i++)
    cFloat[i] = aFloat[i] - bFloat[i];
  return *(zgfGaugeMatrix *)cFloat;
}

double operator~(zgfGaugeMatrix a)
{
  double c = 0.0;
  c += pow(a.c11.re, 2) + pow(a.c12.re, 2) + pow(a.c13.re, 2) + pow(a.c11.im, 2) + pow(a.c12.im, 2) + pow(a.c13.im, 2);
  c += pow(a.c21.re, 2) + pow(a.c22.re, 2) + pow(a.c23.re, 2) + pow(a.c21.im, 2) + pow(a.c22.im, 2) + pow(a.c23.im, 2);
  c += pow(a.c31.re, 2) + pow(a.c32.re, 2) + pow(a.c33.re, 2) + pow(a.c31.im, 2) + pow(a.c32.im, 2) + pow(a.c33.im, 2);
  return c;
}

void zgfGenAField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf)
{
#pragma omp parallel for
  for (int i = 0; i < VOL * Nd; i++)
  {
    double tr = (gf[i].c11.im + gf[i].c22.im + gf[i].c33.im) / 3;
    af[i].c11.im = 0.0;
    af[i].c11.re = gf[i].c11.im - tr;
    af[i].c22.im = 0.0;
    af[i].c22.re = gf[i].c22.im - tr;
    af[i].c33.im = 0.0;
    af[i].c33.re = gf[i].c33.im - tr;

    af[i].c12.im = (gf[i].c12.re - gf[i].c21.re) / -2;
    af[i].c12.re = (gf[i].c12.im + gf[i].c21.im) / 2;
    af[i].c13.im = (gf[i].c13.re - gf[i].c31.re) / -2;
    af[i].c13.re = (gf[i].c13.im + gf[i].c31.im) / 2;
    af[i].c23.im = (gf[i].c23.re - gf[i].c32.re) / -2;
    af[i].c23.re = (gf[i].c23.im + gf[i].c32.im) / 2;

    af[i].c21.re = af[i].c12.re;
    af[i].c21.im = -af[i].c12.im;
    af[i].c31.re = af[i].c13.re;
    af[i].c31.im = -af[i].c13.im;
    af[i].c32.re = af[i].c23.re;
    af[i].c32.im = -af[i].c23.im;
  }
}

void zgfGenDeltaField(zgfGaugeMatrix *df, zgfGaugeMatrix *af)
{
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
  {
    int t = i / (Nz * Ny * Nz);
    int z = i % (Nz * Ny * Nx) / (Ny * Nx);
    int y = i % (Ny * Nx) / (Nx);
    int x = i % (Nx);
    int jt = i + (t == 0 ? (Nt - 1) : -1) * Nz * Ny * Nx;
    int jz = i + (z == 0 ? (Nz - 1) : -1) * Ny * Nx;
    int jy = i + (y == 0 ? (Ny - 1) : -1) * Nx;
    int jx = i + (x == 0 ? (Nx - 1) : -1);

    df[i] = af[i * Nd] + af[i * Nd + 1] + af[i * Nd + 2] + af[i * Nd + 3] - af[jt * 4 + 3] - af[jz * Nd + 2] - af[jy * Nd + 1] - af[jx * Nd];
  }
}

void zgfGenKField(zgfGaugeMatrix *kf, zgfGaugeMatrix *gf, zgfGaugeMatrix *grf)
{
}

double zgfGetTheta(zgfGaugeMatrix *df)
{
  double *theta;
  theta = (double *)malloc(VOL * sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < VOL; i++)
    theta[i] = ~(df[i]);
  for (int i = 1; i < VOL; i++)
    theta[0] += theta[i];
  return theta[0] / Nc / VOL;
}