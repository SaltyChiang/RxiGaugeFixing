/*! \file
 *  \brief Generate a Lambda field for Rxi gauge fixing
 */

#include "../include/chromabase.h"
#include "../include/genlambda.h"

namespace Chroma
{

//! Generate a Lambda field for Rxi gauge fixing
/*!
 * \ingroup gfix
 *
 * Generate a Lambda field, whose Lambda^a ~ i.i.d. N(0, xi) 
 * and sum_x(Lambda(x))=0. The field will be used in R_xi
 * gauge fixing iterations.
 *
 * \param lambda     Lambda field used in Rxi gauge fixing ( Modify )
 * \param xi         Rxi gauge fixing parameter ( Read )
 */

void genLambda(LatticeColorMatrix &lambda,
               const double xi)
{
  bool doGen = true;
  const int Ngen = Nc * Nc - 1;
  multi1d<LatticeReal> random(Ngen);
  multi1d<Real> average(Ngen);

  Complex iReal(0), iImag(0);
  iReal.elem().elem().elem().real() = 1.;
  iImag.elem().elem().elem().imag() = 1.;
  LatticeComplex latReal(iReal), latImag(iImag);

  multi1d<LatticeColorMatrix> gellmann(Ngen);

  if (Ngen == 8)
  {
    gellmann[0] = 0;
    pokeColor(gellmann[0], latReal, 0, 1);
    pokeColor(gellmann[0], latReal, 1, 0);

    gellmann[1] = 0;
    pokeColor(gellmann[1], -latImag, 0, 1);
    pokeColor(gellmann[1], latImag, 1, 0);

    gellmann[2] = 0;
    pokeColor(gellmann[2], latReal, 0, 0);
    pokeColor(gellmann[2], -latReal, 1, 1);

    gellmann[3] = 0;
    pokeColor(gellmann[3], latReal, 0, 2);
    pokeColor(gellmann[3], latReal, 2, 0);

    gellmann[4] = 0;
    pokeColor(gellmann[4], -latImag, 0, 2);
    pokeColor(gellmann[4], latImag, 2, 0);

    gellmann[5] = 0;
    pokeColor(gellmann[5], latReal, 1, 2);
    pokeColor(gellmann[5], latReal, 2, 1);

    gellmann[6] = 0;
    pokeColor(gellmann[6], -latImag, 1, 2);
    pokeColor(gellmann[6], latImag, 2, 1);

    gellmann[7] = 0;
    pokeColor(gellmann[7], latReal, 0, 0);
    pokeColor(gellmann[7], latReal, 1, 1);
    pokeColor(gellmann[7], -2. * latReal, 1, 1);
    gellmann[7] /= sqrt(3.);
  }
  else
  {
    doGen = false;
    QDPIO::cout << "genLambda() for Nc!=3 not implemented. Lambda will be set to 0." << std::endl;
  }

  lambda = 0;
  if (doGen)
    for (int a = 0; a < Ngen; ++a)
    {
      gaussian(random[a]);
      random[a] *= sqrt(xi);
      average[a] = sum(random[a]) / double(Layout::vol());
      random[a] -= average[a];
      random[a] *= sqrt(double(Layout::vol()) / double(Layout::vol() - 1));
      lambda += random[a] * gellmann[a];
    }
}

} // end namespace Chroma
