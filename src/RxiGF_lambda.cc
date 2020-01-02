#include <random>

#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_lambda.h"

#define PI 3.14159265358979323846

namespace RxiGF
{

ColorMatrix Gell_Mann[8];

void genLambdaField(ColorMatrix *lf, double xi)
{
    double testLambda[8];
    for (int i = 0; i < 8; i++)
        testLambda[i] = 0.;

    Gell_Mann[0] = ColorMatrix::Zero();
    Gell_Mann[0](0, 1) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[0](1, 0) = Eigen::dcomplex(0.5, 0.);

    Gell_Mann[1] = ColorMatrix::Zero();
    Gell_Mann[1](0, 1) = Eigen::dcomplex(0., -0.5);
    Gell_Mann[1](1, 0) = Eigen::dcomplex(0., 0.5);

    Gell_Mann[2] = ColorMatrix::Zero();
    Gell_Mann[2](0, 0) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[2](1, 1) = Eigen::dcomplex(-0.5, 0.);

    Gell_Mann[3] = ColorMatrix::Zero();
    Gell_Mann[3](0, 2) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[3](2, 0) = Eigen::dcomplex(0.5, 0.);

    Gell_Mann[4] = ColorMatrix::Zero();
    Gell_Mann[4](0, 2) = Eigen::dcomplex(0., -0.5);
    Gell_Mann[4](2, 0) = Eigen::dcomplex(0., 0.5);

    Gell_Mann[5] = ColorMatrix::Zero();
    Gell_Mann[5](1, 2) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[5](2, 1) = Eigen::dcomplex(0.5, 0.);

    Gell_Mann[6] = ColorMatrix::Zero();
    Gell_Mann[6](1, 2) = Eigen::dcomplex(0., -0.5);
    Gell_Mann[6](2, 1) = Eigen::dcomplex(0., 0.5);

    Gell_Mann[7] = ColorMatrix::Zero();
    Gell_Mann[7](0, 0) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[7](1, 1) = Eigen::dcomplex(0.5, 0.);
    Gell_Mann[7](2, 2) = Eigen::dcomplex(-1., 0.);
    Gell_Mann[7] /= sqrt(3);

    std::mt19937_64 mt_rand(7);
    double rand[(Nc * Nc) / 2 * 2];
    double gaus[(Nc * Nc) / 2 * 2];
    for (int i = 0; i < VOL; i++)
    {
        lf[i] = ColorMatrix::Zero();
        for (int j = 0; j < (Nc * Nc) / 2 * 2; j++)
        {
            rand[j] = mt_rand() / double(UINT64_MAX);
            if (j % 2 == 1)
            {
                gaus[j - 1] = cos(2 * PI * rand[j - 1]) * sqrt(-2 * xi * log(1 - rand[j]));
                gaus[j - 0] = sin(2 * PI * rand[j - 1]) * sqrt(-2 * xi * log(1 - rand[j]));
                testLambda[j - 1] += gaus[j - 1];
                testLambda[j - 0] += gaus[j - 0];
                if (j < (Nc * Nc - 1))
                {
                    lf[i] += gaus[j - 1] * Gell_Mann[j - 1];
                    lf[i] += gaus[j - 0] * Gell_Mann[j - 0];
                }
            }
        }
    }
    for (int i = 0; i < 8; i++)
        printf("%le ", testLambda[i]);
    printf("\n");
}

} // namespace RxiGF