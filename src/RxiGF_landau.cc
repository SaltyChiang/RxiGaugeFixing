#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_landau.h"
#include "../include/RxiGF_linalg_eigen.h"
#include "../include/RxiGF_relax.h"
#include "../include/RxiGF_steepest.h"

namespace RxiGF
{

int LandauGaugeRelax(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax, bool overrelax, double overrelaxParam)
{
  int iterCount = 0;
  double res = 1.0;
  double funcOld, funcNew;

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
        RelaxGaugeRotateField_eigen(grf, gf, tgf, su2_index, cb, overrelax, overrelaxParam);
      } /* end su2_index loop */
    }   /* end cb loop */

    /* Reunitarize */
    ReunitGaugeRotateField_eigen(grf);

    UpdateGaugeField_eigen(tgf, gf, grf);

    funcNew = GetFunctional_eigen(tgf);
    /* Normalized convergence criterion: */
    res = fabs((funcNew - funcOld) / funcNew);
    funcOld = funcNew;

    // std::cout << funcNew << " " << funcOld << std::endl;
    std::cout << "COULGAUGE: iter= " << iterCount
              << "  tgfold= " << funcOld
              << "  tgfnew= " << funcNew
              << "  convar= " << res
              << std::endl;
  } /* end while loop */

  UpdateGaugeField_eigen(gf, grf);

  return iterCount;
}

int LandauGaugeSteepest(ColorMatrix *gf, ColorMatrix *tgf, ColorMatrix *grf, double iterAccu, int iterMax)
{
  int iterCount = 0;
  double theta;

  theta = GetTheta_eigen(deltaField, aField, gf);
  InitGaugeRotateField_eigen(grf);

  while ((theta > iterAccu) && iterCount < iterMax)
  {
    iterCount += 1;

    SteepestGaugeRotateField_eigen(grf, deltaField);

    /* Reunitarize */
    ReunitGaugeRotateField_eigen(grf);

    UpdateGaugeField_eigen(tgf, gf, grf);
    printf("%le\n", theta);
    theta = GetTheta_eigen(deltaField, aField, tgf);
  } /* end while loop */

  UpdateGaugeField_eigen(gf, grf);

  return iterCount;
}

} // namespace RxiGF