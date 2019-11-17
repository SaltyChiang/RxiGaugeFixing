#include <cstdio>

#include "include/zetaGF_time.h"
#include "include/zetaGF_lattice.h"
#include "include/zetaGF_io.h"
#include "include/zetaGF_linalg.h"

zgfGaugeField gaugeField;
zgfGaugeField aField;
zgfGaugeMatrix deltaField[Nt][Nz][Ny][Nx];

int main(int argc, char *argv[])
{
  zgfReadConf(&(gaugeField[0][0][0][0][0]), new char[15]{"data/conf.data"});

  StartTime(1);
  zgfGenAField(&(aField[0][0][0][0][0]), &(gaugeField[0][0][0][0][0]));
  zgfGenDeltaField(&(deltaField[0][0][0][0]), &(aField[0][0][0][0][0]));
  StopTime(1);
  StartTime(2);
  double theta = zgfGetTheta(&(deltaField[0][0][0][0]));
  StopTime(2);
  printf("%le\n", theta);
  PrintTime(1, "DeltaField");
  PrintTime(2, "Theta");
  // printf("%le\n", zgfGaugeMatrixTrace(aField[0][0][0][0][0]) - zgfGaugeMatrixTrace(aField[Nt - 1][0][0][0][0]));
  // printf("%le\n", zgfGaugeMatrixTrace(deltaField[0][0][0][0][0]));
  // printf("%le\n", aField[0][0][0][0][0].c11.re - aField[Nt - 1][0][0][0][0].c11.re);
  // printf("%le\n", deltaField[0][0][0][0][0].c11.re);
  // zgfWriteConf(&(gaugeField[0][0][0][0][0]), new char[14]{"data/conf.out"});
  return 1;
}