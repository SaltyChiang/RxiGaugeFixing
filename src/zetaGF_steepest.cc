#include "../include/zetaGF_lattice.h"
#include "../include/zetaGF_relax.h"
#include "../include/zetaGF_linalg_eigen.h"

namespace ZetaGF
{

void SteepestGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *df)
{
#pragma omp parallel for
    for (int i = 0; i < VOL; i++)
    {
        grf[i] = -df[i];
    }
}

} // namespace ZetaGF