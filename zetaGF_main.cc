#include <cstdio>
#include <cstdlib>

#include "include/zetaGF_lattice.h"
#include "include/zetaGF_io.h"
#include "include/zetaGF_linalg.h"
#include "include/zetaGF_linalg_eigen.h"
#include "include/zetaGF_landau.h"

using namespace ZetaGF;

ColorMatrix *gaugeField;       // gf
ColorMatrix *tempGaugeField;   // tgf
ColorMatrix *aField;           // af
ColorMatrix *deltaField;       // df
ColorMatrix *gaugeRotateField; //grf
// ColorMatrix *kField; // kf
ColorMatrix *v;
double *realA;
double *realB;
double *r_l;
bool *lbtmp;

int main(int argc, char *argv[])
{
  gaugeField = zgfMalloc(ColorMatrix, VOL * Nd);
  tempGaugeField = zgfMalloc(ColorMatrix, VOL * Nd);
  aField = zgfMalloc(ColorMatrix, VOL * Nd);
  deltaField = zgfMalloc(ColorMatrix, VOL);
  gaugeRotateField = zgfMalloc(ColorMatrix, VOL);
  // kField = zgfMalloc(ColorMatrix, VOL);
  v = zgfMalloc(ColorMatrix, VOL * Nd);
  realB = zgfMalloc(double, VOL * 4);
  realA = zgfMalloc(double, VOL * 4);
  r_l = zgfMalloc(double, VOL);
  lbtmp = zgfMalloc(bool, VOL);
  ReadConf(gaugeField, new char[15]{"data/conf.data"});

  StartTimeChrono(1);
  double theta1 = GetTheta_eigen(deltaField, aField, gaugeField);
  StopTimeChrono(1);
  printf("%le\n", theta1);

  LandauGauge(gaugeField, tempGaugeField, gaugeRotateField, 1e-5, 1, true, 1.7);

  double theta2 = GetTheta_eigen(deltaField, aField, gaugeField);
  printf("%le\n", theta2);

  // StartTimeChrono(2);
  // double theta2 = 0.0;
  // double  theta2 = GetTheta_eigen(deltaField, aField, gaugeField);
  // StopTimeChrono(2);

  PrintTimeChrono(1, "Mine");
  // PrintTimeChrono(2, "Eigen");
  return 1;
}