#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_relax.h"
#include "../include/RxiGF_linalg_eigen.h"

namespace RxiGF
{

void SteepestGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *df)
{
#pragma omp parallel for
    for (int i = 0; i < VOL; i++)
    {
        grf[i] = -df[i];
    }
}

} // namespace RxiGF