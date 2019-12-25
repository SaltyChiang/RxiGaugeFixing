#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "extension/Eigen/Dense"

int main()
{
    Eigen::Matrix3cd a[4];
    Eigen::Matrix3cd *b = (Eigen::Matrix3cd *)malloc(4 * 9 * 2 * sizeof(double));
    double *c = (double *)b;

    FILE *fp = fopen("data/test.double", "rb");
    fread((void *)c, sizeof(double) * 2, 4 * 9, fp);
    for (int ilat = 0; ilat < 4; ilat++)
        for (int ic1 = 0; ic1 < 3; ic1++)
            for (int ic2 = 0; ic2 < 3; ic2++)
            {
                a[ilat](ic1, ic2) = Eigen::dcomplex(c[ilat * 18 + ic1 * 6 + ic2 * 2 + 0], c[ilat * 18 + ic1 * 6 + ic2 * 2 + 1]);
            }
    std::cout << a[3] << std::endl;
    std::cout << b[3] << std::endl;

    return 1;
}