#pragma once

#include "RxiGF_lattice.h"

namespace RxiGF
{

void CopyGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf);

void InitGaugeRotateField_eigen(ColorMatrix *grf);
void ReunitGaugeRotateField_eigen(ColorMatrix *grf);

void UpdateGaugeField_eigen(ColorMatrix *gf, ColorMatrix *grf);
void UpdateGaugeField_eigen(ColorMatrix *tgf, ColorMatrix *gf, ColorMatrix *grf);

void GenAField_eigen(ColorMatrix *af, ColorMatrix *gf);
void GenDeltaField_eigen(ColorMatrix *df, ColorMatrix *af);
double GetTheta_eigen(ColorMatrix *df, ColorMatrix *af, ColorMatrix *gf);
double GetTheta_eigen(ColorMatrix *df, ColorMatrix *af, ColorMatrix *gf, ColorMatrix *lf);

double GetFunctional_eigen(ColorMatrix *gf);
double GetFunctional_eigen(ColorMatrix *gf, ColorMatrix *grf, ColorMatrix *lf);

/*
 * Deprecated.
 */
void GenKField_eigen(ColorMatrix *af, ColorMatrix *gf, ColorMatrix *grf);

} // namespace RxiGaugeFixing