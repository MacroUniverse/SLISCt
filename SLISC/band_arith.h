#pragma once
#include "cmat.h"

namespace slisc {

template <class T, SLS_IF(is_scalar<T>())>
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

template <class T, SLS_IF(is_scalar<T>())>
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

template <class T, SLS_IF(is_scalar<T>())>
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

template <class T, SLS_IF(is_scalar<T>())>
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

} // namespace slisc
