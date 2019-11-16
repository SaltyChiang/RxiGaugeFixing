#pragma once

#include "zetaGF_lattice.h"

double zgfGaugeMatrixTrace(zgfGaugeMatrix mat);
void zgfGenAField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf);
void zgfGenDeltaField(zgfGaugeMatrix *df, zgfGaugeMatrix *af);
double zgfGetTheta(zgfGaugeMatrix *df);