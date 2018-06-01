// nr3.h extension
// written by Hongyu Shi
#pragma once
#include "nr3.h"

const Doub PI{ 3.14159265358979323 };
const Complex I(0, 1);
const Doub E{exp(1.)};

// display vectors and matrices
// version 1
void disp(VecDoub_I &v);
void disp(VecComplex_I &v);
void disp(MatDoub_I &v);
void disp(MatComplex_I &v);
// version 2
void disp(VecDoub_I &v, const int precision);
void disp(VecComplex_I &v, const int precision);
void disp(MatDoub_I &v, const int precision);
void disp(MatComplex_I &v, const int precision);
// version 3
void disp(VecDoub_I &v, const int start, const int n);
void disp(VecComplex_I &v, const int start, const int n);
void disp(MatDoub_I &v, const int start1, const int start2, const int n1, const int n2);
void disp(MatComplex_I &v, const int start1, const int start2, const int n1, const int n2);
// version 4
void disp(VecDoub_I &v, const int start, const int n, const int precision);
void disp(VecComplex_I &v, const int start, const int n, const int precision);
void disp(MatDoub_I &v, const int start1, const int start2, const int n1, const int n2, const int precision);
void disp(MatComplex_I &v, const int start1, const int start2, const int n1, const int n2, const int precision);


// numel function returns total number of elements
inline Int numel(VecDoub_I &v) { return v.size(); }
inline Int numel(MatDoub_I &v) { return v.nrows()*v.ncols(); }
inline size_t numel(Mat3DDoub_I &v) { return v.dim1()*(size_t)v.dim2()*v.dim3(); }
inline Int numel(VecComplex_I &v) { return v.size(); }
inline Int numel(MatComplex_I &v) { return v.nrows()*v.ncols(); }
inline size_t numel(Mat3DComplex_I &v) { return v.dim1()*(size_t)v.dim2()*v.dim3(); }
inline Int numel(VecBool_I &v) { return v.size(); }
inline Int numel(MatBool_I &v) { return v.nrows()*v.ncols(); }

// pointer function returns a pointer to the first element
inline const Doub* pointer(VecDoub_I &v) { return &v[0]; }
inline Doub* pointer(VecDoub &v) { return &v[0]; }
inline const Doub* pointer(MatDoub_I &v) { return &v[0][0]; }
inline Doub* pointer(MatDoub &v) { return &v[0][0]; }
inline const Complex* pointer(VecComplex_I &v) { return &v[0]; }
inline Complex* pointer(VecComplex &v) { return &v[0]; }
inline const Complex* pointer(MatComplex_I &v) { return &v[0][0]; }
inline Complex* pointer(MatComplex &v) { return &v[0][0]; }
inline const Complex* pointer(Mat3DComplex_I &v) { return &v[0][0][0]; }
inline Complex* pointer(Mat3DComplex &v) { return &v[0][0][0]; }
inline const Doub* pointer(Mat3DDoub_I &v) { return &v[0][0][0]; }
inline Doub* pointer(Mat3DDoub &v) { return &v[0][0][0]; }

// Matrix/Vector manipulation

template <class T>
inline void flip(NRvector<T> &v)
{
	Int i, n{ v.size() }, ind;
	T temp;
	for (i = 0; i < n / 2; ++i) {
		ind = n - i - 1;
		temp = v[i]; v[i] = v[ind]; v[ind] = temp;
	}
}

template <class T>
inline void flip(NRvector<T> &v, const NRvector<T> &v0)
{
	Int i, n{ v0.size() };
	T temp;
	if (v.size() != n) v.resize(n);
	for (i = 0; i < n; ++i)
		v[i] = v0[n - i - 1];
}

// operators between Complex and Int

inline Complex operator+(const Complex c, const Int i) { return c + (Doub)i; }

inline Complex operator+(const Int i, const Complex c) { return c + (Doub)i; }

inline Complex operator-(const Int i, const Complex c) { return (Doub)i - c; }

inline Complex operator-(const Complex c, const Int i) { return c - (Doub)i; }

inline Complex operator*(const Complex c, const Int i) { return c*(Doub)i; }

inline Complex operator*(const Int i, const Complex c) { return c*(Doub)i; }

inline Complex operator/(const Complex c, const Int i) { return c/(Doub)i; }

inline Complex operator/(const Int i, const Complex c) { return (Doub)i/c; }

// element-wise operators for vectors and matrices
template <class T, class S>
inline void operator<<(NRmatrix<T> &v, const S s)
{
	int i, N{ numel(v) };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] = s;
}

template <class T, class S>
inline void operator<<(NRvector<T> &v, const S s)
{
	int i, N{ numel(v) };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] = s;
}

// +,-,*,/ between vectors and scalars
// for performence, use add(v, v1, v2), minus(v, v1, v2) etc. instead

template <class T>
inline NRvector<T> operator+(const NRvector<T> &v1, const Doub s)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] + s;
}

template <class T>
inline NRvector<T> operator+(const Doub s, const NRvector<T> &v1)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] + s;
}

template <class T>
inline NRvector<T> operator-(const NRvector<T> &v1, const Doub s)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

template <class T>
inline NRvector<T> operator-(const Doub s, const NRvector<T> &v1)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

template <class T>
inline NRvector<T> operator*(const NRvector<T> &v1, const Doub s)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

template <class T>
inline NRvector<T> operator*(const Doub s, const NRvector<T> &v1)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

template <class T>
inline NRvector<T> operator/(const NRvector<T> &v1, const Doub s)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] / s;
}

template <class T>
inline NRvector<T> operator/(const Doub s, const NRvector<T> &v1)
{
	int i, N{ v1.size() };
	NRvector<T> v(N);
	for (i = 0; i < N; ++i)
		v[i] = s / v1[i];
}

template <class T>
inline void operator+=(NRvector<T> &v1, const NRvector<T> & v2)
{
	int i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] += v2[i];
}

template <class T>
inline void operator+=(NRmatrix<T> &v1, const NRmatrix<T> & v2)
{
	int i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] += pv2[i];
}

template <class T>
inline void operator+=(NRMat3d<T> &v1, const NRMat3d<T> & v2)
{
	size_t i, N{ v1.dim1()*(size_t)v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] += pv2[i];
}

template <class T>
inline void operator-=(NRvector<T> &v1, const NRvector<T> & v2)
{
	int i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] -= v2[i];
}

template <class T>
inline void operator-=(NRmatrix<T> &v1, const NRmatrix<T> & v2)
{
	int i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] -= pv2[i];
}

template <class T>
inline void operator-=(NRMat3d<T> &v1, const NRMat3d<T> & v2)
{
	size_t i, N{ v1.dim1()*(size_t)v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] -= pv2[i];
}

template <class T>
inline void operator*=(NRvector<T> &v1, const NRvector<T> & v2)
{
	int i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] *= v2[i];
}

template <class T>
inline void operator*=(NRmatrix<T> &v1, const NRmatrix<T> & v2)
{
	int i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] *= pv2[i];
}

template <class T>
inline void operator*=(NRMat3d<T> &v1, const NRMat3d<T> & v2)
{
	size_t i, N{ v1.dim1()*(size_t)v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] *= pv2[i];
}

template <class T>
inline void operator/=(NRvector<T> &v1, const NRvector<T> & v2)
{
	int i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] /= v2[i];
}

template <class T>
inline void operator/=(NRmatrix<T> &v1, const NRmatrix<T> & v2)
{
	int i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] /= pv2[i];
}

template <class T>
inline void operator/=(NRMat3d<T> &v1, const NRMat3d<T> & v2)
{
	size_t i, N{ v1.dim1()*(size_t)v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] /= pv2[i];
}

template <class T>
inline void add(NRvector<T> &v, const NRvector<T> &v1, const NRvector<T> &v2)
{
	Int i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

template <class T>
inline void add(NRmatrix<T> &v, const NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Int i, N{ v1.nrows()*v1.ncols() };
	auto pv = v[0];
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] + pv2[i];
}

template <class T>
inline void add(NRMat3d<T> &v, const NRMat3d<T> &v1, const NRMat3d<T> & v2)
{
	size_t i, N{ v1.dim1()*(size_t)v1.dim2()*v1.dim3() };
	auto pv = v[0][0];
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] + pv2[i];
}

template <class T>
inline void emul(T &v, const T &v1, const T &v2)
{
	Int i, N{ numel(v1) };
	auto pv = pointer(v);
	auto pv1 = pointer(v1);
	auto pv2 = pointer(v2);
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] * pv2[i];
}

// vector/matrix scalar operator
template <class T>
inline void operator+=(NRvector<T> &v, const Doub s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] += s;
	}
}

template <class T>
inline void operator+=(NRmatrix<T> &v, const Doub s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] += s;
	}
}

template <class T>
inline void operator+=(NRMat3d<T> &v, const Doub s)
{
	size_t i, N{ v.dim1()*(size_t)v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] += s;
	}
}

template <class T>
inline void operator-=(NRvector<T> &v, const Doub s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] -= s;
	}
}

template <class T>
inline void operator-=(NRmatrix<T> &v, const Doub s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] -= s;
	}
}

template <class T>
inline void operator-=(NRMat3d<T> &v, const Doub s)
{
	size_t i, N{ v.dim1()*(size_t)v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] -= s;
	}
}

template <class T>
inline void operator*=(NRvector<T> &v, const Doub s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] *= s;
	}
}

template <class T>
inline void operator*=(NRmatrix<T> &v, const Doub s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] *= s;
	}
}

template <class T>
inline void operator*=(NRMat3d<T> &v, const Doub s)
{
	size_t i, N{ v.dim1()*(size_t)v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] *= s;
	}
}

template <class T>
inline void operator/=(NRvector<T> &v, const Doub s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] /= s;
	}
}

template <class T>
inline void operator/=(NRmatrix<T> &v, const Doub s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] /= s;
	}
}

template <class T>
inline void operator/=(NRMat3d<T> &v, const Doub s)
{
	size_t i, N{ v.dim1()*(size_t)v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] /= s;
	}
}

template <class T>
inline void operator+=(NRvector<T> &v, const Complex s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] += s;
	}
}

template <class T>
inline void operator+=(NRmatrix<T> &v, const Complex s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] += s;
	}
}

template <class T>
inline void operator-=(NRvector<T> &v, const Complex s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] -= s;
	}
}

template <class T>
inline void operator-=(NRmatrix<T> &v, const Complex s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] -= s;
	}
}

template <class T>
inline void operator*=(NRvector<T> &v, const Complex s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] *= s;
	}
}

template <class T>
inline void operator*=(NRmatrix<T> &v, const Complex s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] *= s;
	}
}

template <class T>
inline void operator/=(NRvector<T> &v, const Complex s)
{
	Int i, N{ v.size() };
	for (i = 0; i < N; ++i) {
		v[i] /= s;
	}
}

template <class T>
inline void operator/=(NRmatrix<T> &v, const Complex s)
{
	Int i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] /= s;
	}
}

// dot products
inline Doub operator*(VecDoub_I &v1, VecDoub_I &v2)
{
	Int i, N{ v1.size() };
	Doub s{};
	for (i = 0; i < N; ++i)
		s += v1[i] * v2[i];
	return s;
}

inline Complex operator*(VecComplex_I &v1, VecComplex_I &v2)
{
	Int i, N{ v1.size() };
	Complex s{};
	for (i = 0; i < N; ++i)
		s += conj(v1[i]) * v2[i];
	return s;
}

// matrix vector multiplication
inline void mul(VecDoub_O &y, MatDoub_I &a, VecDoub_I &x)
{
	Int i, j, Nx{ a.nrows() }, Ny{ a.ncols() };
	for (i = 0; i < Ny; ++i) {
		for (j = 0; j < Nx; ++j) {
			y[i] = a[i][j] * x[j];
		}
	}
}

inline void mul(VecComplex_O &y, MatDoub_I &a, VecComplex_I &x)
{
	Int i, j, Nx{ a.ncols() }, Ny{ a.nrows() };
	for (i = 0; i < Ny; ++i) {
		for (j = 0; j < Nx; ++j) {
			y[i] = a[i][j] * x[j];
		}
	}
}

inline void mul(VecComplex_O &y, MatComplex_I &a, VecDoub_I &x)
{
	Int i, j, Nx{ a.ncols() }, Ny{ a.nrows() };
	for (i = 0; i < Ny; ++i) {
		for (j = 0; j < Nx; ++j) {
			y[i] = a[i][j] * x[j];
		}
	}
}

// matrix matrix multiplication
inline void mul(MatDoub_O &c, MatDoub_I &a, MatDoub_I &b)
{
	Int i, j, k, m{ c.nrows() }, n{ c.ncols() }, Nk{ a.ncols() };
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			for (k = 0; k < Nk; ++k) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

// norm and norm square of a vector
inline Doub norm(VecDoub_I &v)
{
	Int i, N{ v.size() };
	Doub s{};
	for (i = 0; i < N; ++i) {
		s += v[i] * v[i];
	}
	return sqrt(s);
}

inline Doub norm2(VecDoub_I &v)
{
	Int i, N{ v.size() };
	Doub s{};
	for (i = 0; i < N; ++i) {
		s += v[i] * v[i];
	}
	return s;
}

inline Doub norm(VecComplex_I &v)
{
	Int i, N{ v.size() };
	Doub s{}, temp;
	for (i = 0; i < N; ++i) {
		temp = abs(v[i]);
		s += temp*temp;
	}
	return sqrt(s);
}

inline Doub norm2(VecComplex_I &v)
{
	Int i, N{ v.size() };
	Doub s{}, temp;
	for (i = 0; i < N; ++i) {
		temp = abs(v[i]);
		s += temp * temp;
	}
	return s;
}

inline Doub norm2(MatComplex_I &v)
{
	Int i, N{ numel(v) };
	Doub s{}, temp;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		temp = abs(pv[i]);
		s += temp * temp;
	}
	return s;
}

inline Doub norm2(Mat3DComplex_I &v)
{
	size_t i, N{ numel(v) };
	Doub s{}, temp;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		temp = abs(pv[i]);
		s += temp * temp;
	}
	return s;
}

// math functions for vectors and matrices

inline void linspace(VecDoub_O &v, const Doub first, const Doub last, Int N = -1)
{
	Int i;
	if (N < 0) N = v.size();
	else v.resize(N);
	Doub delta{ (last - first) / (N - 1) };
	for (i = 0; i < N; ++i) {
		v[i] = first + delta * i;
	}
}

inline void linspace(VecComplex_O &v, const Complex first, const Complex last, Int N = -1)
{
	Int i;
	if (N < 0) N = v.size();
	else v.resize(N);
	Complex delta{ (last - first) / (Doub)(N - 1) };
	for (i = 0; i < N; ++i) {
		v[i] = first + delta * (Doub)i;
	}
}

inline void linspace(MatDoub_O &a, const Doub first, const Doub last)
{
	int i, N{ a.nrows()*a.ncols() };
	Doub delta{ (last - first) / (N - 1) };
	Doub *pa{ &a[0][0] };
	for (i = 0; i < N; ++i) {
		pa[i] = first + delta *i;
	}
}

template <class T>
void exp(NRvector<T> &y, const NRvector<T> &x)
{
	Int i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = exp(x[i]);
	}
}

template <class T>
void exp(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Int i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = exp(px[i]);
	}
}

template <class T>
void sin(NRvector<T> &y, const NRvector<T> &x)
{
	Int i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = sin(x[i]);
	}
}

template <class T>
void sin(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Int i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = sin(px[i]);
	}
}

template <class T>
void cos(NRvector<T> &y, const NRvector<T> &x)
{
	Int i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = cos(x[i]);
	}
}

template <class T>
void cos(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Int i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = cos(px[i]);
	}
}

template <class T>
void tan(NRvector<T> &y, const NRvector<T> &x)
{
	Int i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = tan(x[i]);
	}
}

template <class T>
void tan(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Int i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = tan(px[i]);
	}
}
