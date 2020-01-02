#include <cstdio>
#include <cstdlib>

#include "include/RxiGF_lattice.h"
#include "include/RxiGF_io.h"
#include "include/RxiGF_linalg.h"
#include "include/RxiGF_linalg_eigen.h"
#include "include/RxiGF_gfix.h"
#include "include/RxiGF_lambda.h"

using namespace RxiGF;

ColorMatrix *gaugeField;       // gf
ColorMatrix *tempGaugeField;   // tgf
ColorMatrix *aField;           // af
ColorMatrix *deltaField;       // df
ColorMatrix *gaugeRotateField; //grf
ColorMatrix *lambdaField;      // lf
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
  lambdaField = zgfMalloc(ColorMatrix, VOL);
  v = zgfMalloc(ColorMatrix, VOL * Nd);
  realB = zgfMalloc(double, VOL * 4);
  realA = zgfMalloc(double, VOL * 4);
  r_l = zgfMalloc(double, VOL);
  lbtmp = zgfMalloc(bool, VOL);
  ReadConf(gaugeField, new char[20]{"data/qio.double"});

  genLambdaField(lambdaField, 0.01);

  double theta1 = GetTheta_eigen(deltaField, aField, gaugeField, lambdaField);
  // StopTimeChrono(1);

  printMatrix(aField[0]);
  printMatrix(lambdaField[0]);

  // StartTimeChrono(1);
  // LandauGaugeRelax(gaugeField, tempGaugeField, gaugeRotateField, 1e-10, 100, true, 1.7);
  RxiGaugeRelax(gaugeField, tempGaugeField, gaugeRotateField, lambdaField, 1e-10, 50000, true, 1.7);
  // // LandauGaugeSteepest(gaugeField, tempGaugeField, gaugeRotateField, 1e-5, 100);
  printf("%le\n", theta1);

  double theta2 = GetTheta_eigen(deltaField, aField, gaugeField, lambdaField);
  printf("%le\n", theta2);
  printMatrix(deltaField[0]);
  printMatrix(lambdaField[0]);

  // PrintTimeChrono(1, "Mine");
  return 1;
}