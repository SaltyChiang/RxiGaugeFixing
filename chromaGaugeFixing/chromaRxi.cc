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
  std::string inputFile, outputFile;
  int lx, ly, lz, lt;

  OptionParser parser(argc, argv);
  inputFile = parser.parseOption('i');
  outputFile = parser.parseOption('o');
  lx = atoi(parser.parseOption('x', "24"));
  ly = atoi(parser.parseOption('y', "24"));
  lz = atoi(parser.parseOption('z', "24"));
  lt = atoi(parser.parseOption('t', "64"));
  const double xi = atof(parser.parseOption('r', "0.0"));
  const double gfAccu = atof(parser.parseOption('a', "1e-10"));
  const double gfMax = atoi(parser.parseOption('m', "1000"));
  const bool orDo = atoi(parser.parseOption('d', "0"));
  const double orPara = atof(parser.parseOption('p', "1.7"));

  int n_gf = 0;
  QDP_initialize(&argc, &argv);
  multi1d<int> nsize(Nd);
  const int nsize_int[] = {lx, ly, lz, lt};
  nsize = nsize_int;
  Layout::setLattSize(nsize);
  Layout::create();
  multi1d<LatticeColorMatrix> u(Nd);
  LatticeColorMatrix lambda;

  Chroma::genLambda(lambda, xi);
  Chroma::readKYU(u, inputFile);
  // Chroma::coulGauge(u, n_gf, Nd, gfAccu, gfMax, orDo, orPara);
  Chroma::rxiGauge(u, lambda, n_gf, Nd, gfAccu, gfMax, orDo, orPara);

  QDP_finalize();

  return 1;
}
