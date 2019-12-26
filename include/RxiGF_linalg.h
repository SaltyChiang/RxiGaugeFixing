#pragma once

#include "RxiGF_lattice.h"

namespace RxiGF
{

void GenAField(GaugeMatrix *af, GaugeMatrix *gf);
void GenDeltaField(GaugeMatrix *df, GaugeMatrix *af);
void GenKField(GaugeMatrix *af, GaugeMatrix *gf, GaugeMatrix *grf);
double GetTheta(GaugeMatrix *df);

} // namespace RxiGaugeFixing