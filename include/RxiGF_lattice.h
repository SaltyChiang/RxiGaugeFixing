#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "../extension/Eigen/Core"

#define Nc 3
#define Nd 4
#define Nt 64
#define Ns 24
#define Nx Ns
#define Ny Ns
#define Nz Ns
#define VOL (Nt * Nz * Ny * Nx)
#define PREC double
#define zgfKappa 0.135
#define zgfBeta 0.1
#define zgfG0 sqrt(2 * Nc / beta)

namespace RxiGF
{

typedef struct
{
  PREC re, im;
} Complex;

typedef struct
{
  Complex c11, c12, c13;
  Complex c21, c22, c23;
  Complex c31, c32, c33;
} GaugeMatrix;

typedef Eigen::Matrix<Eigen::dcomplex, Nc, Nc, Eigen::RowMajor> ColorMatrix;

inline void printMatrix(ColorMatrix matrix)
{
  for (int i = 0; i < Nc; i++)
  {
    for (int j = 0; j < Nc; j++)
      printf("(%10.6f,%10.6f) ", matrix(i, j).real(), matrix(i, j).imag());
    printf("\n");
  }
}

} // namespace RxiGF

#include <chrono>
#include "RxiGF_macro.h"

extern RxiGF::ColorMatrix *v;
extern double *realA;
extern double *realB;
extern double *r_l;
extern bool *lbtmp;

extern RxiGF::ColorMatrix *deltaField;
extern RxiGF::ColorMatrix *aField;