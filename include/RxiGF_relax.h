#pragma once

#include "RxiGF_lattice.h"

namespace RxiGF
{

void RelaxGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *gf, ColorMatrix *tgf, int su2_index, int cb, bool overrelax, double overrelaxParam);
void RelaxGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *lf, int su2_index, int cb, bool overrelax, double overrelaxParam);

} // namespace RxiGF