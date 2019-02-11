// display vectors and matrices
// header-only version, in namespace slisc
// unable to be called in debugger

#pragma once
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "fixsize.h"
#include "sparse.h"

namespace slisc {

const Int def_disp_prec = 4;

// ptr version
template <class T>
void disp(const T *ptr, Long_I N, Int_I precision = def_disp_prec);

// version 1 & 2
template <class T>
void disp(const Vector<T> &v, Int_I precision = def_disp_prec);
template <class T>
void disp(const Matrix<T> &a, Int_I precision = def_disp_prec);
template <class T>
void disp(const Cmat<T> &a, Int_I precision = def_disp_prec);
template <class T, Long Nr, Long Nc>
void disp(const FixCmat<T, Nr, Nc> &a, Int_I precision = def_disp_prec);
template <class T>
void disp(const MatCoo<T> &a, Int_I precision = def_disp_prec);
template <class T>
void disp(const MatCooH<T> &a, Int_I precision = def_disp_prec);
template <class T>
void disp(const Mat3d<T> &a, Int_I precision = def_disp_prec);

// version 3 & 4
template <class T>
void disp(const Vector<T> &v, Long_I start, Long_I n, Int_I precision = def_disp_prec);
template <class T>
void disp(const Matrix<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision = def_disp_prec);
template <class T>
void disp(const Cmat<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision = def_disp_prec);
template <class T>
void disp(const MatCoo<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision = def_disp_prec);
template <class T>
void disp(const MatCooH<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision = def_disp_prec);
template <class T>
void disp(const Mat3d<T> &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision = def_disp_prec);

// implementation

// Char2Int<T>::type is Int when T is Char
// otherwise return the same T

// version 0

template <class T>
void disp(const T *ptr, Long_I n, Int_I precision)
{
	for (Int i = 0; i < n; ++i) {
		cout << to_num(ptr[i]) << "   ";
	}
}

// version 1 & 2

template <class T>
inline void disp(const Vector<T> &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
		for (i = 0; i < n; ++i) {
			cout << to_num(v[i]) << "   ";
		}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

template <class T>
inline void disp_mat(const T &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
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

template <class T>
void disp(const Matrix<T> &a, Int_I precision)
{ disp_mat(a, precision); }

template <class T>
void disp(const Cmat<T> &a, Int_I precision)
{ disp_mat(a, precision); }

template <class T, Long Nr, Long Nc>
void disp(const FixCmat<T, Nr, Nc> &a, Int_I precision)
{ disp_mat(a, precision); }

template <class T>
void disp(const MatCoo<T> &a, Int_I precision)
{ disp_mat(a, precision); }

template <class T>
void disp(const MatCooH<T> &a, Int_I precision)
{ disp_mat(a, precision); }

template <class T>
inline void disp(const Mat3d<T> &a, Int_I precision)
{
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

// version 3 & 4

template <class T>
inline void disp(const Vector<T> &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << to_num(v[i]) << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

template <class T>
inline void disp_mat(const T &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << to_num(a(i, j)) << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

template <class T>
inline void disp(const Matrix<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp_mat(a, start1, start2, n1, n2, precision); }

template <class T>
inline void disp(const Cmat<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp_mat(a, start1, start2, n1, n2, precision); }

template <class T>
inline void disp(const MatCoo<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp_mat(a, start1, start2, n1, n2, precision); }

template <class T>
inline void disp(const MatCooH<T> &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp_mat(a, start1, start2, n1, n2, precision); }

template <class T>
inline void disp(const Mat3d<T> &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{
	Long i, j, k;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (k = start3; k < start3 + n3; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				cout << a(i, j, k) << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

} // namespace slisc
