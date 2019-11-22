#pragma once

#include "zetaGF_lattice.h"

namespace ZetaGF
{

void RelaxGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *gf, ColorMatrix *tgf, int su2_index, int cb, bool overrelax, double overrelaxParam);

} // namespace ZetaGF