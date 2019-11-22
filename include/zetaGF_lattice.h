#pragma once

#include "../extension/Eigen/Core"

#define Nc 3
#define Nd 4
#define Nt 32
#define Ns 32
#define Nx Ns
#define Ny Ns
#define Nz Ns
#define VOL (Nt * Nz * Ny * Nx)
#define PREC double
#define zgfKappa 0.135
#define zgfBeta 0.1
#define zgfG0 sqrt(2 * Nc / beta)

namespace ZetaGF
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

typedef Eigen::Matrix<Eigen::dcomplex, Nc, Nc> ColorMatrix;

} // namespace ZetaGF

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include "zetaGF_macro.h"

extern ZetaGF::ColorMatrix *v;
extern double *realA;
extern double *realB;
extern double *r_l;
extern bool *lbtmp;

extern ZetaGF::ColorMatrix *deltaField;
extern ZetaGF::ColorMatrix *aField;
