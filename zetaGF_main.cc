#include <cstdio>
#include <cstdlib>

#include "include/zetaGF_lattice.h"
#include "include/zetaGF_io.h"
#include "include/zetaGF_linalg.h"
#include "include/zetaGF_linalg_eigen.h"
#include "include/zetaGF_landau.h"

using namespace ZetaGF;

ColorMatrix *gaugeField; // gf
ColorMatrix *tempGaugeField; // tgf
ColorMatrix *aField; // af
ColorMatrix *deltaField; // df
ColorMatrix *gaugeRotateField; //grf
ColorMatrix *kField; // kf

int main(int argc, char *argv[])
{
  gaugeField = zgfMalloc(ColorMatrix, VOL * Nd);
  tempGaugeField = zgfMalloc(ColorMatrix, VOL * Nd);
  aField = zgfMalloc(ColorMatrix, VOL * Nd);
  deltaField = zgfMalloc(ColorMatrix, VOL);
  gaugeRotateField = zgfMalloc(ColorMatrix, VOL);
  kField = zgfMalloc(ColorMatrix, VOL);
  ReadConf(gaugeField, new char[15]{"data/conf.data"});

  StartTimeChrono(1);
  // zgfGenAField_eigen(aField, gaugeField);
  // zgfGenDeltaField_eigen(deltaField, aField);
  // double theta1 = zgfGetTheta_eigen(deltaField);
  // zgfInitGaugeRotateField_eigen(gaugeRotateField);
  StopTimeChrono(1);
  // printf("%le\n", theta1);

  // StartTimeChrono(2);
  // double theta2 = 0.0;
  // for (int i = 0; i < 5; i++)
  // {
  //   zgfGenKField_eigen(kField, gaugeField, gaugeRotateField);
  //   zgfUpdateGaugeField_eigen(gaugeField, gaugeRotateField);
  //   zgfGenAField_eigen(aField, gaugeField);
  //   zgfGenDeltaField_eigen(deltaField, aField);
  //   theta2 = zgfGetTheta_eigen(deltaField);
  //   printf("%le\n", theta2);
  // }
  // StopTimeChrono(2);

  PrintTimeChrono(1, "Mine");
  // PrintTimeChrono(2, "Eigen");
  LandauGauge(gaugeField, tempGaugeField, gaugeRotateField, 1e-10, 100, false, 1.7);
  return 1;
}