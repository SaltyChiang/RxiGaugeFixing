#include <cstdio>
#include <cstdlib>

#include "include/chromabase.h"
#include "include/genlambda.h"
#include "include/coulgauge.h"
#include "include/rxigauge.h"
#include "include/gauge_io.h"
#include "include/kyugauge_io.h"
#include "include/helpfunc.h"

int main(int argc, char *argv[])
{
  OptionParser parser(argc, argv);

  int n_gf = 0;
  QDP_initialize(&argc, &argv);
  multi1d<int> nsize(Nd);
  const int nsize_int[] = {parser.x, parser.y, parser.z, parser.t};
  nsize = nsize_int;
  Layout::setLattSize(nsize);
  Layout::create();
  multi1d<LatticeColorMatrix> u(Nd);
  LatticeColorMatrix lambda;

  QDPIO::cout << "xi = " << parser.xi << std::endl
              << "GFAccu = " << parser.gfAccu << std::endl
              << "GFMax = " << parser.gfMax << std::endl
              << "OrDo = " << parser.orDo << std::endl
              << "OrPara = " << parser.orPara << std::endl;

  Chroma::genLambda(lambda, parser.xi);
  Chroma::readKYU(u, parser.input);
  // Chroma::coulGauge(u, n_gf, Nd, gfAccu, gfMax, orDo, orPara);
  Chroma::rxiGauge(u, lambda, n_gf, Nd, parser.gfAccu, parser.gfMax, parser.orDo, parser.orPara);

  QDP_finalize();

  return 1;
}
