#pragma once

#include "zetaGF_lattice.h"

namespace ZetaGF
{

void CopyGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf);
void GenAField_eigen(ColorMatrix *af, ColorMatrix *gf);
void GenDeltaField_eigen(ColorMatrix *df, ColorMatrix *af);
void InitGaugeRotateField_eigen(ColorMatrix *grf);
void ReunitGaugeRotateField_eigen(ColorMatrix *grf);
void GenKField_eigen(ColorMatrix *af, ColorMatrix *gf, ColorMatrix *grf);
void GenGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *tgf, int su2_index, int cb, bool overrelax, double overrelaxParam);
void UpdateGaugeField_eigen(ColorMatrix *gf, ColorMatrix *grf);
void UpdateGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf, ColorMatrix *grf);
double GetTheta_eigen(ColorMatrix *df);
double GetFunctional_eigen(ColorMatrix *gf);

} // namespace ZetaGaugeFixing