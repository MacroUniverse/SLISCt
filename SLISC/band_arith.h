#pragma once
#include "band.h"

namespace slisc {

// conversion between full matrix and band diagonal matrix
// Nup and Nlow are the number of upper diagonals and lower diagonals
// ref: cBLAS gbmv() routine
// https://software.intel.com/en-us/node/834918#DAEC7CD0-620A-4696-9612-C295F8211646
template <class T, SLS_IF0(is_scalar<T>())>
void mat2band(Cmat<T> &b, const Cmat<T> &a, Long_I Nup, Long_I Nlow)
{
    Long N1 = a.n1(), N2 = a.n2();
    for (Long j = 0; j < N2; j++) {
        Long k = Nup - j;
        for (Long i = MAX(0, j - Nup); i < MIN(N1, j + Nlow + 1); i++) {
            b(k + i, j) = a(i, j);
        }
    }
}

template <class T, SLS_IF0(is_scalar<T>())>
void mat2band(Matrix<T> &b, const Matrix<T> &a, Long_I Nup, Long_I Nlow)
{
    Long N1 = a.n1(), N2 = a.n2();
    for (Long i = 0; i < N1; i++) {
        Long k = Nlow - i;
        for (Long j = MAX(0, i - Nlow); j < MIN(N2, i + Nup + 1); j++) {
            b(i, k + j) = a(i, j);
        }
    }
}

template <class T, SLS_IF0(is_scalar<T>())>
void band2mat(Cmat<T> &a, const Cmat<T> &b, Long_I Nup, Long_I Nlow)
{
    Long N1 = a.n1(), N2 = a.n2();
    for (Long j = 0; j < N2; j++) {
        Long k = Nup - j;
        for (Long i = MAX(0, j - Nup); i < MIN(N1, j + Nlow + 1); i++) {
            a(i, j) = b(k + i, j);
        }
    }
}

template <class T, SLS_IF0(is_scalar<T>())>
void band2mat(Matrix<T> &a, const Matrix<T> &b, Long_I Nup, Long_I Nlow)
{
    Long N1 = a.n1(), N2 = a.n2();
    for (Long i = 0; i < N1; i++) {
        Long k = Nlow - i;
        for (Long j = MAX(0, i - Nlow); j < MIN(N2, i + Nup + 1); j++) {
            a(i, j) = b(i, k + j);
        }
    }
}

// matrix-vector multiplication for band matrix

// matrix-vector multiplication for band matrix
void mul_gb(VecDoub_O y, const Band<Doub> &a, VecDoub_I x)
{
    Doub alpha = 1, beta = 0;
    Long lda = a.nlow() + a.nup() + 1;
    Long incx = 1, incy = 1;
    cblas_dgbmv(CblasColMajor, CblasNoTrans, a.n1(), a.n2(), a.nlow(), a.nup(),
        alpha, a.ptr(), lda, x.ptr(), incx, beta, y.ptr(), incy);
}

void mul_gb(VecComp_O y, const Band<Comp> &a, VecComp_I x)
{
    Comp alpha(1, 0), beta(0, 0);
    Long lda = a.nlow() + a.nup() + 1;
    Long incx = 1, incy = 1;
    cblas_zgbmv(CblasColMajor, CblasNoTrans, a.n1(), a.n2(), a.nlow(), a.nup(),
        &alpha, a.ptr(), lda, x.ptr(), incx, &beta, y.ptr(), incy);
}

} // namespace slisc
