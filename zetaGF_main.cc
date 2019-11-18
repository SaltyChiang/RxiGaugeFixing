#include <cstdio>
#include <cstdlib>

#include "include/zetaGF_time.h"
#include "include/zetaGF_time_chrono.h"
#include "include/zetaGF_lattice.h"
#include "include/zetaGF_io.h"
#include "include/zetaGF_linalg.h"
#include "include/zetaGF_linalg_eigen.h"

zgfGaugeMatrix *gaugeField;
zgfGaugeMatrix *aField;
zgfGaugeMatrix *deltaField;
zgfGaugeMatrix *gaugeRotateField;
zgfGaugeMatrix *kField;

int main(int argc, char *argv[])
{
  gaugeField = (zgfGaugeMatrix *)malloc(sizeof(zgfGaugeMatrix) * VOL * Nd);
  aField = (zgfGaugeMatrix *)malloc(sizeof(zgfGaugeMatrix) * VOL * Nd);
  deltaField = (zgfGaugeMatrix *)malloc(sizeof(zgfGaugeMatrix) * VOL);
  gaugeRotateField = (zgfGaugeMatrix *)malloc(sizeof(zgfGaugeMatrix) * VOL);
  kField = (zgfGaugeMatrix *)malloc(sizeof(zgfGaugeMatrix) * VOL);
  zgfReadConf(gaugeField, new char[15]{"data/conf.data"});
  printf("???\n");

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