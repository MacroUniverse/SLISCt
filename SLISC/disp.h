// display vectors and matrices
// header-only version, in namespace slisc
// unable to be called in debugger

#pragma once
#include "scalar_arith.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "cmat3d.h"
#include "fixsize.h"
#include "matcooh.h"
#include "cmatobd.h"

namespace slisc {

const Int def_disp_prec = 8;

// display vector
template <class Tv, SLS_IF(ndims<Tv>() == 1)>
void disp(const Tv &v, Int_I precision = def_disp_prec)
{
    auto oldPrecision = cout.precision();
    cout.precision(precision);
    if (v.size() == 0)
        cout << "empty";
    else
        for (Long i = 0; i < v.size(); ++i) {
            cout << to_num(v[i]) << "   ";
        }
    cout << endl << endl;
    cout.precision(oldPrecision);
}

// display matrix
template <class T, SLS_IF(ndims<T>() == 2)>
void disp(const T &a, Int_I precision = def_disp_prec)
{
    Long i, j, m{ a.n1() }, n{ a.n2() };
    auto oldPrecision = cout.precision();
    cout.precision(precision);
    if (a.size() == 0) cout << "empty";
    else
        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
                cout << to_num(a(i, j)) << "   ";
            }
            cout << endl;
        }
    cout << endl;
    cout.precision(oldPrecision);
}

// display 3d array
template <class T, SLS_IF(ndims<T>() == 3)>
void disp(const T &a, Int_I precision = def_disp_prec) {
    Long i, j, k, m{ a.n1() }, n{ a.n2() }, q{ a.n3() };
    auto oldPrecision = cout.precision();
    cout.precision(precision);
    if (a.size() == 0) cout << "empty";
    else
        for (k = 0; k < q; ++k) {
            cout << "(:, :, " << k << ")" << endl;
            for (i = 0; i < m; ++i) {
                for (j = 0; j < n; ++j) {
                    cout << to_num(a(i, j, k)) << "   ";
                }
                cout << endl;
            }
            cout << endl;
        }
    cout.precision(oldPrecision);
}

// display 4d array
template <class T, SLS_IF(ndims<T>() == 4)>
void disp(const T &a, Int_I precision = def_disp_prec) {
    Long i, j, k, l, N1 = a.n1(), N2 = a.n2(), N3 = a.n3(), N4 = a.n4();
    auto oldPrecision = cout.precision();
    cout.precision(precision);
    if (a.size() == 0)
        cout << "empty";
    else
        for (l = 0; l < N4; ++l)
        for (k = 0; k < N3; ++k) {
            cout << "(:, :, " << k << ", " << l << ")" << endl;
            for (i = 0; i < N1; ++i) {
                for (j = 0; j < N2; ++j) {
                    cout << to_num(a(i, j, k, l)) << "   ";
                }
                cout << endl;
            }
            cout << endl;
        }
    cout.precision(oldPrecision);
}
} // namespace slisc
