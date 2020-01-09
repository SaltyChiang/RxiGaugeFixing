#pragma once

#include "RxiGF_lattice.h"

namespace RxiGF
{

int LandauGaugeRelax(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam);
int RxiGaugeRelax(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, ColorMatrix *lf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam);

int LandauGaugeSteepest(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax);

} // namespace RxiGF