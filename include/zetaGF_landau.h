#pragma once

#include "zetaGF_lattice.h"

namespace ZetaGF
{

int LandauGaugeRelax(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam);

int LandauGaugeSteepest(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax);

} // namespace ZetaGF