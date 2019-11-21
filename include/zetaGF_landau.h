#pragma once

#include "zetaGF_lattice.h"

namespace ZetaGF
{

int LandauGauge(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam);

} // namespace ZetaGaugeFixing