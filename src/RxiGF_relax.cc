#include "../include/RxiGF_lattice.h"
#include "../include/RxiGF_relax.h"
#include "../include/RxiGF_linalg_eigen.h"

namespace RxiGF
{

int i1, i2;
int found;
int del_i;
int index;
const double fuzz = 1e-10;

void su2Extract(double *r, ColorMatrix *c, int su2_index)
{
    /* Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k */
    r[0] = (*c)(i1, i1).real() + (*c)(i2, i2).real();
    r[1] = (*c)(i1, i2).imag() + (*c)(i2, i1).imag();
    r[2] = (*c)(i1, i2).real() - (*c)(i2, i1).real();
    r[3] = (*c)(i1, i1).imag() - (*c)(i2, i2).imag();
}

void sunFill(double *r, ColorMatrix *c, int su2_index)
{

    /*
     * Insert the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k
     * back into the SU(N) matrix
     */
    *c = ColorMatrix::Identity();

    (*c)(i1, i1) = Eigen::dcomplex(r[0], r[3]);
    (*c)(i1, i2) = Eigen::dcomplex(r[2], r[1]);
    (*c)(i2, i1) = Eigen::dcomplex(-r[2], r[1]);
    (*c)(i2, i2) = Eigen::dcomplex(r[0], -r[3]);
}

void gRelax(ColorMatrix *tgf, int su2_index, int cb, int i, int xm, int ym, int zm, int tm, bool overrelax, double overrelaxParam)
{
    v[i] = tgf[i * Nd + 0] + tgf[xm * Nd + 0].adjoint();
    v[i] += tgf[i * Nd + 1] + tgf[ym * Nd + 1].adjoint();
    v[i] += tgf[i * Nd + 2] + tgf[zm * Nd + 2].adjoint();
    v[i] += tgf[i * Nd + 3] + tgf[tm * Nd + 3].adjoint();

    su2Extract(&(realB[i * Nd]), &(v[i]), su2_index);

    /*
     * Now project onto SU(2)
     */
    r_l[i] = sqrt(realB[i * Nd + 0] * realB[i * Nd + 0] + realB[i * Nd + 1] * realB[i * Nd + 1] + realB[i * Nd + 2] * realB[i * Nd + 2] + realB[i * Nd + 3] * realB[i * Nd + 3]);

    // Normalize
    lbtmp[i] = r_l[i] > 1e-10;

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    realA[i * Nd + 0] = lbtmp[i] ? (realB[i * Nd + 0] / r_l[i]) : 1.0;
    realA[i * Nd + 1] = lbtmp[i] ? -(realB[i * Nd + 1] / r_l[i]) : 0.0;
    realA[i * Nd + 2] = lbtmp[i] ? -(realB[i * Nd + 2] / r_l[i]) : 0.0;
    realA[i * Nd + 3] = lbtmp[i] ? -(realB[i * Nd + 3] / r_l[i]) : 0.0;

    /* Now do the overrelaxation, if desired */
    if (overrelax)
    {
        /* get angle */
        double theta_old;
        theta_old = acos(realA[i * Nd + 0]);

        /* old sin */
        double oldsin;
        oldsin = sin(theta_old);

        /* overrelax, i.e. multiply by the angle */
        double theta_new;
        theta_new = theta_old * overrelaxParam;

        /* compute sin(new)/sin(old) */
        /* set the ratio to 0, if sin(old) < FUZZ */
        //lftmp = (oldsin > fuzz) ? (sin(theta_new) / oldsin) : 0.0;

        /* get the new cos = a[0] */
        realA[i * Nd + 0] = cos(theta_new);

        /* get the new a_k, k = 1, 2, 3 */
        realA[i * Nd + 1] *= (oldsin > fuzz) ? (sin(theta_new) / oldsin) : 0.0;
        realA[i * Nd + 2] *= (oldsin > fuzz) ? (sin(theta_new) / oldsin) : 0.0;
        realA[i * Nd + 3] *= (oldsin > fuzz) ? (sin(theta_new) / oldsin) : 0.0;
    }

    /* Now fill the SU(Nc) matrix V with the SU(2) submatrix 'su2_index' */
    /* paramtrized by a_k in the sigma matrix basis. */
    sunFill(&(realA[i * 4]), &(v[i]), su2_index);
}

void RelaxGaugeRotateField_eigen(ColorMatrix *grf, ColorMatrix *gf, ColorMatrix *tgf, int su2_index, int cb, bool overrelax, double overrelaxParam)
{
    UpdateGaugeField_eigen(tgf, gf, grf);

    found = 0;
    del_i = 0;
    index = -1;

    while ((del_i < (Nc - 1)) && (found == 0))
    {
        del_i++;
        for (i1 = 0; i1 < (Nc - del_i); i1++)
        {
            index++;
            if (index == su2_index)
            {
                found = 1;
                break;
            }
        }
    }
    i2 = i1 + del_i;

    if (found == 0)
    {
        exit(255);
    }

#pragma omp parallel for
    for (int i = 0; i < VOL; i++)
    {
        v[i] = ColorMatrix::Zero();
        int t = i / (Nz * Ny * Nz);
        int z = i % (Nz * Ny * Nx) / (Ny * Nx);
        int y = i % (Ny * Nx) / (Nx);
        int x = i % (Nx);
        int tm = i + (t == 0 ? (Nt - 1) : -1) * Nz * Ny * Nx;
        int zm = i + (z == 0 ? (Nz - 1) : -1) * Ny * Nx;
        int ym = i + (y == 0 ? (Ny - 1) : -1) * Nx;
        int xm = i + (x == 0 ? (Nx - 1) : -1);
        // int tp = i + (t == (Nt - 1) ? (1 - Nt) : 1) * Nz * Ny * Nx;
        // int zp = i + (z == (Nz - 1) ? (1 - Nz) : 1) * Ny * Nx;
        // int yp = i + (y == (Ny - 1) ? (1 - Ny) : 1) * Nx;
        // int xp = i + (x == (Nx - 1) ? (1 - Nx) : 1);
        if ((t + z + y + x) % 2 == cb)
        {
            gRelax(tgf, su2_index, cb, i, xm, ym, zm, tm, overrelax, overrelaxParam);
            grf[i] = v[i] * grf[i];
        }
    }
}

} // namespace RxiGF