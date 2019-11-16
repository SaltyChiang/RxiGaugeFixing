#include <cstdlib>

#include "../include/zetaGF_lattice.h"
#include "../include/zetaGF_linalg.h"

zgfGaugeMatrix zgfGaugeMatrixSub(zgfGaugeMatrix a, zgfGaugeMatrix b)
{
  PREC *cFloat = (PREC *)malloc(sizeof(zgfGaugeMatrix));
  PREC *aFloat = (PREC *)&a;
  PREC *bFloat = (PREC *)&b;
  for (int i = 0; i < Nc * Nc * 2; i++)
    cFloat[i] = aFloat[i] - bFloat[i];
  return *(zgfGaugeMatrix *)cFloat;
}

double zgfGaugeMatrixTrace(zgfGaugeMatrix mat)
{
  return mat.c11.re + mat.c22.re + mat.c33.re;
}

void zgfGenAField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf)
{
#pragma omp parallel for
  for (int i = 0; i < V * Nd; i++)
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
  for (int i = 0; i < V * Nd; i++)
  {
    int t = i / (Nz * Ny * Nz * Nd);
    int z = i % (Nz * Ny * Nx * Nd) / (Ny * Nx * Nd);
    int y = i % (Ny * Nx * Nd) / (Nx * Nd);
    int x = i % (Nx * Nd) / (Nd);
    int d = i % (Nd);
    int a = 0;
    int j = i;
    switch (d)
    {
    case 0:
      a = (t == 0 ? Nt - 1 : t - 1);
      j += (a - t) * Nz * Ny * Nx * Nd;
      break;
    case 1:
      a = (z == 0 ? Nz - 1 : z - 1);
      j += (a - z) * Ny * Nx * Nd;
      break;
    case 2:
      a = (y == 0 ? Ny - 1 : y - 1);
      j += (a - y) * Nx * Nd;
      break;
    case 3:
      a = (x == 0 ? Nx - 1 : x - 1);
      j += (a - x) * Nd;
      break;
    }
    df[i] = zgfGaugeMatrixSub(af[i], af[j]);
  }
}

double zgfGetTheta(zgfGaugeMatrix *df)
{
  return 1.0;
}