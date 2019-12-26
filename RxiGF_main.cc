#include <cstdio>
#include <cstdlib>

#include "include/RxiGF_lattice.h"
#include "include/RxiGF_io.h"
#include "include/RxiGF_linalg.h"
#include "include/RxiGF_linalg_eigen.h"
#include "include/RxiGF_landau.h"

using namespace RxiGF;

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
  ReadConf(gaugeField, new char[20]{"data/qio.double"});
  // printMatrix(gaugeField[0]);
  // printMatrix(gaugeField[VOL * Nd - 1]);

  // StartTimeChrono(1);
  double theta1 = GetTheta_eigen(deltaField, aField, gaugeField);
  // StopTimeChrono(1);
  printf("%le\n", theta1);

  LandauGaugeRelax(gaugeField, tempGaugeField, gaugeRotateField, 1e-10, 1000, true, 1.7);
  // // LandauGaugeSteepest(gaugeField, tempGaugeField, gaugeRotateField, 1e-5, 100);

  double theta2 = GetTheta_eigen(deltaField, aField, gaugeField);
  printf("%le\n", theta2);

  // PrintTimeChrono(1, "Mine");
  return 1;
}