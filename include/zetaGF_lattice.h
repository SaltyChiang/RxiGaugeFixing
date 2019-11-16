#pragma once

#define Nc 3
#define Nd 4
#define Nt 32
#define Ns 32
#define Nx Ns
#define Ny Ns
#define Nz Ns
#define V Nt *Nz *Ny *Nx
#define PREC double
#define kappa 0.135
#define beta 0.1
#define g0 sqrt(2 * Nc / beta)

typedef struct
{
  PREC re, im;
} zgfComplex;

typedef struct
{
  zgfComplex c11, c12, c13;
  zgfComplex c21, c22, c23;
  zgfComplex c31, c32, c33;
} zgfGaugeMatrix;

typedef zgfGaugeMatrix zgfGaugeField[Nt][Nz][Ny][Nx][Nd];

extern zgfGaugeField deltaField;