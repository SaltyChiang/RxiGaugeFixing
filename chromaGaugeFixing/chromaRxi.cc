#include <cstdio>
#include <cstdlib>

#include "qdp.h"
#include "qdp_iogauge.h"

using namespace QDP;
#include "include/coulgauge.h"
#include "include/commio.h"

const int precision = sizeof(double);

int ReadConf(multi1d<LatticeColorMatrix> &gf, char *fn)
{
  FILE *fp = fopen(fn, "rb");
  double *ptr = (double *)malloc(Layout::vol() * Nd * Nc * Nc * 2 * precision);
  fread((void *)ptr, 2 * precision, Layout::vol() * Nd * Nc * Nc, fp);
  // printf("%le\n", ptr[0]);

  for (int i = 0; i < Layout::vol(); i++)
  {
    int t = i / (16 * 16 * 16);
    int z = i % (16 * 16 * 16) / (16 * 16);
    int y = i % (16 * 16) / (16);
    int x = i % (16);
    int cb = (x + y + z + t) % 2;
    x = x / 2;
    int ii = t * (16 * 16 * 8) + z * (16 * 8) + y * 8 + x;
    for (int id = 0; id < Nd; id++)
      for (int ic1 = 0; ic1 < Nc; ic1++)
        for (int ic2 = 0; ic2 < Nc; ic2++)
        {
          gf[id].elem(cb * 32768 + ii).elem().elem(ic1, ic2).real() = ptr[i * Nd * Nc * Nc * 2 + id * Nc * Nc * 2 + ic1 * Nc * 2 + ic2 * 2 + 0];
          gf[id].elem(cb * 32768 + ii).elem().elem(ic1, ic2).imag() = ptr[i * Nd * Nc * Nc * 2 + id * Nc * Nc * 2 + ic1 * Nc * 2 + ic2 * 2 + 1];
        }
  }
  fclose(fp);
  free(ptr);
  return 1;
}

int main(int argc, char *argv[])
{
  int n_gf = 0;

  QDP_initialize(&argc, &argv);
  multi1d<int> nsize(Nd);
  const int nsize_int[] = {16, 16, 16, 16};
  nsize = nsize_int;
  Layout::setLattSize(nsize);
  Layout::create();
  multi1d<LatticeColorMatrix> u(Nd);

  ReadConf(u, (char *)"../data/qio.double");
  // for (int i = 0; i < Nd; i++)
  //   gaussian(u[i]);
  // printMatrix(u[0].elem(0).elem());
  // printMatrix(u[Nd - 1].elem(Layout::vol() - 1).elem());

  Chroma::coulGauge(u, n_gf, Nd, 1e-10, 1000, true, 1.7);

  QDP_finalize();

  return 1;
}