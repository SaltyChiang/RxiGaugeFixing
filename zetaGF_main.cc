#include <cstdio>

#include "include/zetaGF_macro.h"
#include "include/zetaGF_lattice.h"
#include "include/zetaGF_io.h"
#include "include/zetaGF_linalg.h"
#include "include/zetaGF_linalg_eigen.h"

zgfColorMatrix *gaugeField;
zgfColorMatrix *aField;
zgfColorMatrix *deltaField;
zgfColorMatrix *gaugeRotateField;
zgfColorMatrix *kField;

int main(int argc, char *argv[])
{
  gaugeField = zgfMalloc(zgfColorMatrix, VOL * Nd);
  aField = zgfMalloc(zgfColorMatrix, VOL * Nd);
  deltaField = zgfMalloc(zgfColorMatrix, VOL);
  gaugeRotateField = zgfMalloc(zgfColorMatrix, VOL);
  kField = zgfMalloc(zgfColorMatrix, VOL);
  zgfReadConf(gaugeField, new char[15]{"data/conf.data"});

  StartTimeChrono(1);
  zgfGenAField_eigen(aField, gaugeField);
  zgfGenDeltaField_eigen(deltaField, aField);
  double theta1 = zgfGetTheta_eigen(deltaField);
  zgfInitGaugeRotateField_eigen(gaugeRotateField);
  StopTimeChrono(1);
  printf("%le\n", theta1);

  StartTimeChrono(2);
  double theta2 = 0.0;
  for (int i = 0; i < 5; i++)
  {
    zgfGenKField_eigen(kField, gaugeField, gaugeRotateField);
    zgfGenGaugeRotateField_eigen(gaugeRotateField, kField);
    zgfUpdateGaugeField_eigen(gaugeField, gaugeRotateField);
    zgfGenAField_eigen(aField, gaugeField);
    zgfGenDeltaField_eigen(deltaField, aField);
    theta2 = zgfGetTheta_eigen(deltaField);
    printf("%le\n", theta2);
  }
  StopTimeChrono(2);

  PrintTimeChrono(1, "Mine");
  PrintTimeChrono(2, "Eigen");
  return 1;
}