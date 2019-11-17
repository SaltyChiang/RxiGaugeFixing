#pragma once

#include "zetaGF_lattice.h"

void zgfGenAField_eigen(zgfGaugeMatrix *af, zgfGaugeMatrix *gf);
void zgfGenDeltaField_eigen(zgfGaugeMatrix *df, zgfGaugeMatrix *af);
void zgfInitGaugeRotateField_eigen(zgfGaugeMatrix *grf);
void zgfGenKField_eigen(zgfGaugeMatrix *af, zgfGaugeMatrix *gf, zgfGaugeMatrix *grf);
void zgfGenGaugeRotateField_eigen(zgfGaugeMatrix *grf, zgfGaugeMatrix *kf);
void zgfUpdateGaugeField_eigen(zgfGaugeMatrix *gf, zgfGaugeMatrix *grf);
double zgfGetTheta_eigen(zgfGaugeMatrix *df);