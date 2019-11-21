#include "../include/zetaGF_lattice.h"
#include "../include/zetaGF_landau.h"
#include "../include/zetaGF_linalg_eigen.h"

namespace ZetaGF
{

int LandauGauge(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam)
{
  int iterCount = 0;
  double res = 1.0;
  double *theta;
  double funcOld, funcNew;
  theta = (double *)malloc(VOL * Nd * sizeof(double));

  funcOld = GetFunctional_eigen(gf);
  InitGaugeRotateField_eigen(grf);

  while ((res > iterAccu) && iterCount < iterMax)
  {
    iterCount += 1;

    /* Loop over checkerboards for gauge fixing */
    for (int cb = 0; cb < 2; ++cb)
    {
      /* Loop over SU(2) subgroup index */
      for (int su2_index = 0; su2_index < Nc * (Nc - 1) / 2; ++su2_index)
      {
        /* Now do a gauge fixing relaxation step */
        GenGaugeRotateField_eigen(grf, tgf, su2_index, cb, overrelax, overrelaxParam);
      } /* end su2_index loop */
    }   /* end cb loop */

    /* Reunitarize */
    ReunitGaugeRotateField_eigen(grf);

    UpdateGaugeField_eigen(tgf, gf, grf);
    funcNew = GetFunctional_eigen(tgf);

    /* Normalized convergence criterion: */
    res = fabs((funcNew - funcOld) / funcNew);
    funcOld = funcNew;
    printf("%le\n", res);
  } /* end while loop */

  UpdateGaugeField_eigen(gf, grf);

  return iterCount;
}

} // namespace ZetaGaugeFixing