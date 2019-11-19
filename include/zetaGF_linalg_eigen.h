#pragma once

#include "zetaGF_lattice.h"

void zgfGenAField_eigen(zgfColorMatrix *af, zgfColorMatrix *gf);
void zgfGenDeltaField_eigen(zgfColorMatrix *df, zgfColorMatrix *af);
void zgfInitGaugeRotateField_eigen(zgfColorMatrix *grf);
void zgfGenKField_eigen(zgfColorMatrix *af, zgfColorMatrix *gf, zgfColorMatrix *grf);
void zgfGenGaugeRotateField_eigen(zgfColorMatrix *grf, zgfColorMatrix *kf);
void zgfUpdateGaugeField_eigen(zgfColorMatrix *gf, zgfColorMatrix *grf);
double zgfGetTheta_eigen(zgfColorMatrix *df);