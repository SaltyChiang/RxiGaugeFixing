#include <random>

#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_lambda.h"
#include "../include/RxiGF_io.h"

#define PI 3.14159265358979323846

namespace RxiGF
{

ColorMatrix Gell_Mann[8];

void genNormalWithSum0(double *result, double xi)
{
    double x, y;
    double sum[8];
    const double factor = sqrt(double(VOL) / double(VOL - 1));
    std::mt19937_64 mt_rand(7);
    for (int i = 0; i < VOL * 4; ++i)
    {
        x = mt_rand() / double(UINT64_MAX);
        y = mt_rand() / double(UINT64_MAX);

        result[2 * i + 0] = cos(2.0 * PI * x) * sqrt(-2.0 * xi * log(y));
        result[2 * i + 1] = sin(2.0 * PI * x) * sqrt(-2.0 * xi * log(y));
    }
    for (int i = 0; i < 8; i++)
    {
        sum[i] = 0;
        for (int j = 0; j < VOL; j++)
            sum[i] += result[j * 8 + i];
        sum[i] /= double(VOL);
        for (int j = 0; j < VOL; j++)
        {
            result[j * 8 + i] -= sum[i];
            result[j * 8 + i] *= factor;
        }
    }
    WriteRandom(result, new char[20]{"data/random.data"});
}

void genLambdaField(ColorMatrix *lf, double xi)
{
    Gell_Mann[0] = ColorMatrix::Zero();
    Gell_Mann[0](0, 1) = Eigen::dcomplex(1., 0.);
    Gell_Mann[0](1, 0) = Eigen::dcomplex(1., 0.);
    // Gell_Mann[0] /= 2;

    Gell_Mann[1] = ColorMatrix::Zero();
    Gell_Mann[1](0, 1) = Eigen::dcomplex(0., -1.);
    Gell_Mann[1](1, 0) = Eigen::dcomplex(0., 1.);
    // Gell_Mann[1] /= 2;

    Gell_Mann[2] = ColorMatrix::Zero();
    Gell_Mann[2](0, 0) = Eigen::dcomplex(1., 0.);
    Gell_Mann[2](1, 1) = Eigen::dcomplex(-1., 0.);
    // Gell_Mann[2] /= 2;

    Gell_Mann[3] = ColorMatrix::Zero();
    Gell_Mann[3](0, 2) = Eigen::dcomplex(1., 0.);
    Gell_Mann[3](2, 0) = Eigen::dcomplex(1., 0.);
    // Gell_Mann[3] /= 2;

    Gell_Mann[4] = ColorMatrix::Zero();
    Gell_Mann[4](0, 2) = Eigen::dcomplex(0., -1.);
    Gell_Mann[4](2, 0) = Eigen::dcomplex(0., 1.);
    // Gell_Mann[4] /= 2;

    Gell_Mann[5] = ColorMatrix::Zero();
    Gell_Mann[5](1, 2) = Eigen::dcomplex(1., 0.);
    Gell_Mann[5](2, 1) = Eigen::dcomplex(1., 0.);
    // Gell_Mann[5] /= 2;

    Gell_Mann[6] = ColorMatrix::Zero();
    Gell_Mann[6](1, 2) = Eigen::dcomplex(0., -1.);
    Gell_Mann[6](2, 1) = Eigen::dcomplex(0., 1.);
    // Gell_Mann[6] /= 2;

    Gell_Mann[7] = ColorMatrix::Zero();
    Gell_Mann[7](0, 0) = Eigen::dcomplex(1., 0.);
    Gell_Mann[7](1, 1) = Eigen::dcomplex(1., 0.);
    Gell_Mann[7](2, 2) = Eigen::dcomplex(-2., 0.);
    Gell_Mann[7] /= sqrt(3);
    // Gell_Mann[7] /= 2;

    double *result = zgfMalloc(double, VOL * 8);
    genNormalWithSum0(result, xi);

    for (int i = 0; i < VOL; i++)
    {
        lf[i] = ColorMatrix::Zero();
        for (int j = 0; j < (Nc * Nc - 1); j++)
            lf[i] += result[i * 8 + j] * Gell_Mann[j];
    }
    free(result);
}

} // namespace RxiGF