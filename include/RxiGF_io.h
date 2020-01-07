#pragma once

#include "RxiGF_lattice.h"

namespace RxiGF
{

int ReadConf(ColorMatrix *gf, char *fn);
int WriteConf(ColorMatrix *gf, char *fn);
int WriteLambda(ColorMatrix *lf, char *fn);
int WriteRandom(double *rd, char *fn);

} // namespace RxiGaugeFixing