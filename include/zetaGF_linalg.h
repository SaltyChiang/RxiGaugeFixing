#pragma once

#include "zetaGF_lattice.h"

void zgfGenAField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf);
void zgfGenDeltaField(zgfGaugeMatrix *df, zgfGaugeMatrix *af);
void zgfGenKField(zgfGaugeMatrix *af, zgfGaugeMatrix *gf, zgfGaugeMatrix *grf);
double zgfGetTheta(zgfGaugeMatrix *df);