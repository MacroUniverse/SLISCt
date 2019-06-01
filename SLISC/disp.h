// display vectors and matrices
// header-only version, in namespace slisc
// unable to be called in debugger

#pragma once
#include "scalar_arith.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "fixsize.h"
#include "sparse.h"

namespace slisc {

const Int def_disp_prec = 4;

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
	Long i, j, m{ a.n1() }, n{ a.ncols() };
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
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
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
} // namespace slisc
