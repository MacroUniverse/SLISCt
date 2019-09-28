#pragma once
#include "global.h"
#include "matrix.h"

namespace slisc {
    // default: shift columns to the right n times (n < 0 shift to left)
// column at the end shifts to the other end
// dim = 2: shift rows down (n < 0 shift up)
template <class T>
void shift(Matrix<T> &a, Llong nshift, Int_I dim = 1)
{
    Long Nr = a.n1(), Nc = a.n2(), n;
    if (dim == 2) {
        // I actually want n to be shift to the left
        if (nshift < 0)
            n = (-nshift) % Nc;
        else if (nshift > 0)
            n = Nc - (nshift%Nc);
        else
            return;
        if (n == 0 || n == Nc) return;
        Long i;
        Long sz = n*sizeof(T), sz_ = (Nc-n)*sizeof(T);
        T *temp = new T[n];
        for (i = 0; i < Nr; ++i) {
            memcpy(temp, a[i], sz);
            memcpy(a[i], a[i] + n, sz_);
            memcpy(a[i] + Nc-n, temp, sz);
        }
        delete[] temp;
    }
    else {
        // I actually want n to be shift up
        if (nshift < 0)
            n = (-nshift) % Nr;
        else if (nshift > 0)
            n = Nr - (nshift%Nr);
        else
            return;
        if (n == 0 || n == Nr) return;
        Long sz = n*Nc*sizeof(T);
        Long sz_ = (Nr-n)*Nc*sizeof(T);
        T *temp = new T[n];
        memcpy(temp, a.ptr(), sz);
        memcpy(a.ptr(), a[n], sz_);
        memcpy(a[Nr-n], temp, sz);
        delete[] temp;
    }
}

// shift the i-th line i times to the left, moving the diagonal to the first column
template <class T>
void diagonals(Matrix<T> &a)
{
    Long i, Nr{ a.n1() }, Nc{ a.n2() };
    T *temp = new T[Nc];
    Long szT = sizeof(T);
    for (i = 1; i < Nr; ++i) {
        memcpy(temp, a[i], i*szT);
        memcpy(a[i], a[i] + i, (Nc-i)*szT);
        memcpy(a[i] + Nc-i, temp, i*szT);
    }
    delete[] temp;
}

// parallel version
template <class T>
void diagonals_par(Matrix<T> &a)
{
    SLS_ERR("TODO");
    /*Long i, Nr{ a.n1() }, Nc{ a.n2() };
    Long szT = sizeof(T);
#pragma omp parallel for
    for (i = 1; i < Nr; ++i) {
        T *temp = new T[Nc];
        memcpy(temp, a[i], i*szT);
        memcpy(a[i], a[i] + i, (Nc - i)*szT);
        memcpy(a[i] + Nc - i, temp, i*szT);
        delete[] temp;
    }*/
}

template <class T>
void idiagonals(Matrix<T> &a)
{
    Long i, Nr{ a.n1() }, Nc{ a.n2() };
    T *temp = new T[Nc];
    Long szT = sizeof(T);
    for (i = 1; i < Nr; ++i) {
        memcpy(temp, a[i], (Nc-i)*szT);
        memcpy(a[i], a[i] + (Nc-i), i*szT);
        memcpy(a[i] + i, temp, (Nc-i)*szT);
    }
    delete[] temp;
}

// parallel version
template <class T>
void idiagonals_par(Matrix<T> &a)
{
    SLS_ERR("TODO");
    /*Long i, Nr{ a.n1() }, Nc{ a.n2() };
    Long szT = sizeof(T);
#pragma omp parallel for
    for (i = 1; i < Nr; ++i) {
        T *temp = new T[Nc];
        memcpy(temp, a[i], (Nc - i)*szT);
        memcpy(a[i], a[i] + (Nc - i), i*szT);
        memcpy(a[i] + i, temp, (Nc - i)*szT);
        delete[] temp;
    }*/
}

} // namespace slisc
