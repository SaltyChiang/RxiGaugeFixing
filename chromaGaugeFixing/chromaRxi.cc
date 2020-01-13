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
  int n_gf = 0;

  OptionParser parser(argc, argv);
  inputFile = parser.parseOption('i');
  outputFile = parser.parseOption('o');
  lx = atoi(parser.parseOption('x'));
  ly = atoi(parser.parseOption('y'));
  lz = atoi(parser.parseOption('z'));
  lt = atoi(parser.parseOption('t'));
  const double xi = atof(parser.parseOption('r'));

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
  // Chroma::coulGauge(u, n_gf, Nd, 1e-10, 1000, false, 1.7);
  Chroma::rxiGauge(u, lambda, n_gf, Nd, 1e-10, 1000, false, 1.7);

  QDP_finalize();

  return 1;
}
