#pragma once

#include "zetaGF_lattice.h"

double zgfGaugeMatrixTrace(zgfGaugeMatrix mat);
void zgfGenAField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf);
void zgfGenDeltaField(zgfGaugeMatrix *df, zgfGaugeMatrix *af);
void zgfGenKField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf, zgfGaugeMatrix *grf);
double zgfGetTheta(zgfGaugeMatrix *df);