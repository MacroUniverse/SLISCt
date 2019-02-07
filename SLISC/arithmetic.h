#pragma once
#include "slisc.h"
#include "fixsize.h"
#include "ptr_arith.h"

namespace slisc {

// === get vec/mat properties ===

// check if vec/mat sizes are the same
template <class T1, class T2, class Enable = typename is_slisc2<T1, T2>::type>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
	if constexpr (T1::ndims() == 1 && T2::ndims() == 1) {
		return v1.size() == v2.size();
	}
	else if constexpr (T1::ndims() == 2 && T2::ndims() == 2) {
		return v1.nrows() == v2.nrows() && v1.ncols() == v2.ncols();
	}
	else if constexpr (T1::ndims() == 3 && T2::ndims() == 3) {
		return v1.dim1() == v2.dim1() && v1.dim2() == v2.dim2()
			&& v1.dim3() == v2.dim3();
	}
	return false;
}

// operator== for slisc containers

template <class T1, class T2, class Enable = typename is_slisc2<T1, T2>::type>
Bool operator==(const T1 &v1, const T2 &v2)
{
	return shape_cmp(v1, v2) && equals_to_vv(v1.ptr(), v2.ptr(), v2.size());
}

template <class Tv, class Ts, class Enable = typename is_slisc_other<Tv, Ts>::type>
Bool operator==(const Tv &v, const Ts &s)
{
	equals_to_vs(v.ptr(), s, v.size());
}

template <class Tv, class Ts, class Enable = typename is_slisc_other<Tv, Ts>::type>
Bool operator==(const Ts &s, const Tv &v)
{ return v == s; }

// operator!= for slisc containers
template <class T1, class T2, class Enable = typename has_1_slisc2<T1, T2>::type>
Bool operator!=(const T1 &v1, const T2 &v2)
{
	return !(v1 == v2);
}

// sum of container elements

TODO: disable for MatCooH

template <class T>
inline const T sum(const Vbase<T> &v)
{
	Long i, n{ v.size() };
	T sum = 0;
	for (i = 0; i < n; ++i)
		sum += v(i);
	return sum;
}

inline Long count(VecBool_I v)
{
	Long i, n{ v.size() };
	Long count = 0;
	for (i = 0; i < n; ++i)
		count += v(i);
	return count;
}

template <class T>
inline const T max(const Vbase<T> &v)
{
	Long i, N{ v.size() };
	T val{ v(0) };
	for (i = 1; i < N; ++i)
		if (v(i) > val)
			val = v(i);
	return val;
}

// return max(abs(a(:))
template <class T>
inline const auto max_abs(const Vbase<T> &v)
{
	typename rm_complex<T>::type s{ abs(v[0]) };
	for (Long i = 1; i < v.size(); ++i)
		if (s < abs(v[i]))
			s = abs(v(i));
	return s;
}

template <class T>
inline const T max(Long_O ind, const Vbase<T> &v)
{
	Long i, N{ v.size() };
	T val{ v(0) };
	for (i = 1; i < N; ++i)
		if (v(i) > val) {
			val = v[i]; ind = i;
		}
	return val;
}

inline Doub max(Long_O ind, const Vbase<Comp> &v)
{
	Long i, N{ v.size() };
	Doub val{ abs(v[0]) };
	for (i = 1; i < N; ++i)
		if (abs(v[i]) > val) {
			val = abs(v[i]); ind = i;
		}
	return val;
}

// sum(v(:).^2) for real numbers
template <class T>
const T norm2(Vbase<T> &v)
{
	Long i, N{ v.size() };
	T s2{};
	for (i = 0; i < N; ++i)
		s2 += SQR(v[i]);
	return s2;
}

template <class T>
const T norm(Vbase<T> &v)
{ return std::sqrt(norm2(v)); }

//sum(abs(v(:)). ^ 2) for complex numbers
inline Doub norm2(Vbase<Comp> &v)
{
	Long i, N{ v.size() };
	Doub s2{};
	for (i = 0; i < N; ++i)
		s2 += SQR(abs(v[i]));
	return s2;
}

inline Doub norm(Vbase<Comp> &v)
{ return std::sqrt(norm2(v)); }

// === matrix manipulation ===

// does not work for integers
template <class T, class T1, class T2>
inline void linspace(Vbase<T> &v, const T1 &first, const T2 &last)
{
	Long i, N{ v.size() };
	T delta = (last - first) / T(N - 1);
	for (i = 0; i < N; ++i)
		v[i] = first + delta * T(i);
}

template <class T, class T1, class T2>
inline void linspace(Vector<T> &v, const T1 &first, const T2 &last, Llong_I n)
{ v.resize(n); linspace(v, first, last); }

template <class T, class T1, class T2>
inline void linspace(Matrix<T> &v, const T1 &first, const T2 &last, Llong_I rows, Llong_I cols)
{ v.resize(rows, cols); linspace(v, first, last); }

template <class T, class T1, class T2>
inline void linspace(Mat3d<T> &v, const T first, const T last, Llong_I dim1, Llong_I dim2, Llong_I dim3)
{ v.resize(dim1, dim2, dim3); linspace(v, first, last); }

// element-wise operators for vectors and matrices

// TODO: transpose

// hermitian conjugate
inline void her(MatComp_O h, MatComp_I a)
{
	Long i, j, m = a.nrows(), n = a.ncols();
	h.resize(n, m);
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			h(j, i) = conj(a(i, j));
}

template <class T>
inline void flip(Vector<T> &v)
{
	Long i, n{ v.size() }, ind;
	T temp;
	for (i = 0; i < n / 2; ++i) {
		ind = n - i - 1;
		temp = v[i]; v[i] = v[ind]; v[ind] = temp;
	}
}

template <class T>
inline void flip(Vector<T> &v, const Vector<T> &v0)
{
	Long i, n{ v0.size() };
	v.resize(n);
	for (i = 0; i < n; ++i)
		v[i] = v0[n - i - 1];
}

// default: shift columns to the right n times (n < 0 shift to left)
// column at the end shifts to the other end
// dim = 2: shift rows down (n < 0 shift up)
template <class T>
void shift(Matrix<T> &a, Llong nshift, Int_I dim = 1)
{
	Long Nr = a.nrows(), Nc = a.ncols(), n;
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
		delete temp;
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
		delete temp;
	}
}

// shift the i-th line i times to the left, moving the diagonal to the first column
template <class T>
void diagonals(Matrix<T> &a)
{
	Long i, Nr{ a.nrows() }, Nc{ a.ncols() };
	T *temp = new T[Nc];
	Long szT = sizeof(T);
	for (i = 1; i < Nr; ++i) {
		memcpy(temp, a[i], i*szT);
		memcpy(a[i], a[i] + i, (Nc-i)*szT);
		memcpy(a[i] + Nc-i, temp, i*szT);
	}
	delete temp;
}

// parallel version
template <class T>
void diagonals_par(Matrix<T> &a)
{
	Long i, Nr{ a.nrows() }, Nc{ a.ncols() };
	Long szT = sizeof(T);
	#pragma omp parallel for
	for (i = 1; i < Nr; ++i) {
		T *temp = new T[Nc];
		memcpy(temp, a[i], i*szT);
		memcpy(a[i], a[i] + i, (Nc-i)*szT);
		memcpy(a[i] + Nc-i, temp, i*szT);
		delete temp;
	}
}

template <class T>
void idiagonals(Matrix<T> &a)
{
	Long i, Nr{ a.nrows() }, Nc{ a.ncols() };
	T *temp = new T[Nc];
	Long szT = sizeof(T);
	for (i = 1; i < Nr; ++i) {
		memcpy(temp, a[i], (Nc-i)*szT);
		memcpy(a[i], a[i] + (Nc-i), i*szT);
		memcpy(a[i] + i, temp, (Nc-i)*szT);
	}
	delete temp;
}

template <class T>
void idiagonals_par(Matrix<T> &a)
{
	Long i, Nr{ a.nrows() }, Nc{ a.ncols() };
	Long szT = sizeof(T);
	#pragma omp parallel for
	for (i = 1; i < Nr; ++i) {
		T *temp = new T[Nc];
		memcpy(temp, a[i], (Nc-i)*szT);
		memcpy(a[i], a[i] + (Nc-i), i*szT);
		memcpy(a[i] + i, temp, (Nc-i)*szT);
		delete temp;
	}
}

// === vectorized math functions ===

template <class T, class T1>
void sqrt(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); sqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void sqrt(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); sqrt_vv(v, v1); }

template <class T, class T1>
void sqrt(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); sqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void invSqrt(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); invSqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void invSqrt(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); invSqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void invSqrt(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); invSqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void sin(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); sin_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void sin(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); sin_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void sin(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); sin_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void cos(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); cos_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void cos(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); cos_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void cos(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); cos_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void exp(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); exp_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void exp(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); exp_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void exp(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); exp_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void tan(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); tan_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void tan(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); tan_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1>
void tan(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); tan_vv(v.ptr(), v1.ptr(), v1.size()); }

// === matrix arithmetics ===

template <class T, class T1>
inline void plus_equals_vv(T &v, const T1 &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	plus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void operator+=(Vector<T> &v, const Vector<T1> &v1)
{ plus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator+=(Matrix<T> &v, const Matrix<T1> &v1)
{ plus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator+=(Cmat<T> &v, const Cmat<T1> &v1)
{ plus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator+=(Mat3d<T> &v, const Mat3d<T1> &v1)
{ plus_equals_vv(v, v1); }

// v -= v

template <class T, class T1>
inline void minus_equals_vv(T &v, const T1 &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	minus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void operator-=(Vector<T> &v, const Vector<T1> &v1)
{ minus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator-=(Matrix<T> &v, const Matrix<T1> &v1)
{ minus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator-=(Cmat<T> &v, const Cmat<T1> &v1)
{ minus_equals_vv(v, v1); }

template <class T, class T1>
inline void operator-=(Mat3d<T> &v, const Mat3d<T1> &v1)
{ minus_equals_vv(v, v1); }

// v *= v

template <class T, class T1>
inline void times_equals_vv(T &v, const T1 &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	times_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void operator*=(Vector<T> &v, const Vector<T1> &v1)
{ times_equals_vv(v, v1); }

template <class T, class T1>
inline void operator*=(Matrix<T> &v, const Matrix<T1> &v1)
{ times_equals_vv(v, v1); }

template <class T, class T1>
inline void operator*=(Cmat<T> &v, const Cmat<T1> &v1)
{ times_equals_vv(v, v1); }

template <class T, class T1>
inline void operator*=(Mat3d<T> &v, const Mat3d<T1> &v1)
{ times_equals_vv(v, v1); }

// v /= v

template <class T, class T1>
inline void devide_equals_vv(T &v, const T1 &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	divide_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void operator/=(Vector<T> &v, const Vector<T1> &v1)
{ devide_equals_vv(v, v1); }

template <class T, class T1>
inline void operator/=(Matrix<T> &v, const Matrix<T1> &v1)
{ devide_equals_vv(v, v1); }

template <class T, class T1>
inline void operator/=(Cmat<T> &v, const Cmat<T1> &v1)
{ devide_equals_vv(v, v1); }

template <class T, class T1>
inline void operator/=(Mat3d<T> &v, const Mat3d<T1> &v1)
{ devide_equals_vv(v, v1); }

// v += s
template <class T, class T1>
inline void operator+=(Vector<T> &v, const T1 &s)
{ plus_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator+=(Matrix<T> &v, const T1 &s)
{ plus_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator+=(Mat3d<T> &v, const T1 &s)
{ plus_equals_vs(v.ptr(), s, v.size()); }

// v -= s

template <class T, class T1>
inline void operator-=(Vector<T> &v, const T1 &s)
{ minus_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator-=(Matrix<T> &v, const T1 &s)
{ minus_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator-=(Mat3d<T> &v, const T1 &s)
{ minus_equals_vs(v.ptr(), s, v.size()); }

// v *= s

template <class T, class T1>
inline void operator*=(Vector<T> &v, const T1 &s)
{ times_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator*=(Matrix<T> &v, const T1 &s)
{ times_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator*=(Cmat<T> &v, const T1 &s)
{ times_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator*=(Mat3d<T> &v, const T1 &s)
{ times_equals_vs(v.ptr(), s, v.size()); }

// v /= s

template <class T, class T1>
inline void operator/=(Vector<T> &v, const T1 &s)
{ divide_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator/=(Matrix<T> &v, const T1 &s)
{ divide_equals_vs(v.ptr(), s, v.size()); }

template <class T, class T1>
inline void operator/=(Mat3d<T> &v, const T1 &s)
{ divide_equals_vs(v.ptr(), s, v.size()); }

// TODO: operator /= for integers

// v %= s
template <class T>
inline void operator%=(Vbase<T> &v, const T &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) %= s;
}

// rem(v, s)

template <class T>
inline void rem(Vector<T> &v, const Vector<T> &v1, const T &s)
{ v.resize(v1); rem(v.ptr(), v1.ptr(), s, v1.size()); }

template <class T>
inline void rem(Matrix<T> &v, const Matrix<T> &v1, const T &s)
{ v.resize(v1); rem(v.ptr(), v1.ptr(), s, v1.size()); }

template <class T>
inline void rem(Mat3d<T> &v, const Mat3d<T> &v1, const T &s)
{ v.resize(v1); rem(v.ptr(), v1.ptr(), s, v1.size()); }

// TODO : rem(v, s, v1)
// TODO : rem(v, v1, v2)

// mod(v, v, s)

template <class T>
inline void mod(Vector<T> &v, const Vector<T> &v1, const T &s)
{ v.resize(v1); mod(v.ptr(), v1.ptr(), s, v1.size()); }

template <class T>
inline void mod(Matrix<T> &v, const Matrix<T> &v1, const T &s)
{ v.resize(v1); mod(v.ptr(), v1.ptr(), s, v1.size()); }

template <class T>
inline void mod(Mat3d<T> &v, const Mat3d<T> &v1, const T &s)
{ v.resize(v1); mod(v.ptr(), v1.ptr(), s, v1.size()); }

// TODO : mod(v, s, v1)
// TODO : mod(v, v1, v2)

// plus(v, v, s)

template <class T, class T1, class T2>
inline void plus_vvs(T &v, const T1 &v1, const T2 &s)
{
	v.resize(v1); plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ plus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ plus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ plus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ plus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ plus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ plus_vvs(v, v1, s); }

// plus(v, v, v)

template <class T, class T1, class T2>
inline void plus_vvv(T &v, const T1 &v1, const T2 &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); plus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{ plus_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{ plus_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ plus_vvv(v, v1, v2); }

// minus(v)
template <class T>
inline void minus(Vbase<T> &v)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] *= -1;
}

// minus(v, v)

template <class T, class T1>
inline void minus_vv(T &v, const T1 &v1)
{
	v.resize(v1); minus_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void minus(Vector<T> &v, const Vector<T1> &v1)
{ minus_vv(v, v1); }

template <class T, class T1>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1)
{ minus_vv(v, v1); }

template <class T, class T1>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1)
{ minus_vv(v, v1); }

// minus(v, s, v)

template <class T, class T1, class T2>
inline void minus_vsv(T &v, const T1 &s, const T2 &v1)
{
	v.resize(v1); minus_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ minus_vsv(v, s, v1); }

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ minus_vsv(v, s, v1); }

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ minus_vsv(v, s, v1); }

// minus(v, v, s)

template <class T, class T1, class T2>
inline void minus_vvs(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{
	v.resize(v1); minus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ minus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ minus_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ minus_vvs(v, v1, s); }

// minus(v, v, v)

template <class T, class T1, class T2>
inline void minus_vvv(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); minus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{ minus_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{ minus_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ minus_vvv(v, v1, v2); }

// times(v, v, s)

template <class T, class T1, class T2>
inline void times_vvs(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{
	v.resize(v1); times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class T2>
inline void times(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ times_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ times_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ times_vvs(v, v1, s); }

// times(v, s, v)
template <class T, class T1, class T2>
inline void times(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ times(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ times(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ times(v, v1, s); }

// times(v, v, v)

template <class T, class T1, class T2>
inline void times_vvv(T &v, const T1 &v1, const T2 &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); times_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

template <class T, class T1, class T2>
inline void times(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{ times_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void times(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{ times_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void times(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ times_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void divide_vvs(T &v, const T1 &v1, const T2 &s)
{ v.resize(v1); divide_vvs(v.ptr(), v1.ptr(), s, v1.size()); }

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ divide_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ divide_vvs(v, v1, s); }

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ divide_vvs(v, v1, s); }

// divide(v, s, v)

template <class T, class T1, class T2>
inline void divide_vsv(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{
	v.resize(v1); divide_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ divide_vsv(v, s, v1); }

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ divide_vsv(v, s, v1); }

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ divide_vsv(v, s, v1); }

// divide(v, v, v)

template <class T, class T1, class T2>
inline void divide_vvv(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); divide_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{ divide_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{ divide_vvv(v, v1, v2); }

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ divide_vvv(v, v1, v2); }

// real(v)
inline void real(Vbase<Comp> &v)
{
	Long i, N{ 2 * v.size() };
	Doub *pd = (Doub *)v.ptr();
	for (i = 1; i < N; i += 2)
		pd[i] = 0.;
}

// resl(v, v)

template <class T>
inline void real_vv(T &v, const Vector<Comp> &v1)
{
	v.resize(v1); real_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T>
inline void real(Vector<T> &v, const Vector<Comp> &v1)
{ real_vv(v, v1); }

template <class T>
inline void real(Matrix<T> &v, const Matrix<Comp> &v1)
{ real_vv(v, v1); }

template <class T>
inline void real(Mat3d<T> &v, const Mat3d<Comp> &v1)
{ real_vv(v, v1); }

// imag(v)
inline void imag(Vbase<Comp> &v)
{
	Long i, N{ 2 * v.size() };
	Doub *pd = (Doub *)v.ptr();
	for (i = 0; i < N; i += 2)
		pd[i] = 0.;
}

// imag(v, v)

template <class T>
inline void imag_vv(T &v, const Vector<Comp> &v1)
{
	v.resize(v1); imag_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T>
inline void imag(Vector<T> &v, const Vector<Comp> &v1)
{ imag_vv(v, v1); }

template <class T>
inline void imag(Matrix<T> &v, const Matrix<Comp> &v1)
{ imag_vv(v, v1); }

template <class T>
inline void imag(Mat3d<T> &v, const Mat3d<Comp> &v1)
{ imag_vv(v, v1); }

// abs(v)
template <class T>
inline void abs(Vbase<T> &v)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; i += 2)
		v[i] = std::abs(v[i]);
}

// abs(v, v)

template <class T, class T1>
inline void abs_vv(T &v, const T1 &v1)
{
	v.resize(v1); abs_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1>
inline void abs(Vector<T> &v, const Vector<T1> &v1)
{ abs_vv(v, v1); }

template <class T, class T1>
inline void abs(Matrix<T> &v, const Matrix<T1> &v1)
{ abs_vv(v, v1); }

template <class T, class T1>
inline void abs(Mat3d<T> &v, const Mat3d<T1> &v1)
{ abs_vv(v, v1); }

// doubl2comp(v, v)

template <class T, class T1>
inline void to_comp_vv(T &v, const T1 &v1)
{
	v.resize(v1);
	to_comp_vv(v.ptr(), v1.ptr(), v1.size());
}

inline void to_comp_vv(Vector<Comp> &v, const Vector<Doub> &v1)
{ to_comp_vv(v, v1); }

inline void to_comp_vv(Matrix<Comp> &v, const Matrix<Doub> &v1)
{ to_comp_vv(v, v1); }

inline void to_comp_vv(Mat3d<Comp> &v, const Mat3d<Doub> &v1)
{ to_comp_vv(v, v1); }

// conj(v)

template <class T>
inline void conj(Vbase<T> &v)
{ return; }

template <class T>
inline void conj(Vbase<std::complex<T>> &v)
{
	Long i, N{ 2 * v.size() };
	Doub *p = (Doub *)v.ptr();
	for (i = 1; i < N; i += 2)
		p[i] = -p[i];
}

// dot products ( sum conj(v1[i])*v2[i] )
// s = dot(v, v)

template <class T1, class T2>
inline auto dot(const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	typename std::common_type<T1, T2>::type s;
	dot_svv(s, v1.ptr(), v2.ptr(), v2.size());
	return s;
}

// outer product ( conj(v1[i})*v2[j] )
// outprod(v, v, v)
template <class T, class T1, class T2>
inline void outprod(Matrix<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
	Long i, j, N1{ v1.size() }, N2{ v2.size() };
	Comp *pc, v1_i;
	v.resize(N1, N2);
	for (i = 0; i < N1; ++i) {
		pc = v[i];
		v1_i = v1[i];
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}
}

// parallel version
template <class T, class T1, class T2>
inline void outprod_par(Matrix<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
	Long i, N1{ v1.size() }, N2{ v2.size() };
	v.resize(N1, N2);
	#pragma omp parallel for
	for (i = 0; i < N1; ++i) {
		Long j;
		Comp *pc, v1_i;
		pc = v[i];
		v1_i = v1[i];
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}
}

template<class T, class T2>
inline void outprod(Matrix<T> &v, VecComp_I v1, const Vector<T2> &v2)
{
	Long i, j, N1{ v1.size() }, N2{ v2.size() };
	Comp *pc, v1_i;
	v.resize(N1, N2);
	for (i = 0; i < N1; ++i) {
		pc = v[i];
		v1_i = conj(v1[i]);
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}
}

template<class T, class T2>
inline void outprod_par(Matrix<T> &v, VecComp_I v1, const Vector<T2> &v2)
{
	Long i, N1{ v1.size() }, N2{ v2.size() };
	v.resize(N1, N2);
	#pragma omp parallel for
	for (i = 0; i < N1; ++i) {
		Long j;
		Comp *pc, v1_i;
		pc = v[i];
		v1_i = conj(v1[i]);
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}
}

template <class T, class T1, class T2>
inline void mul_vmv(T &y, const T1 &a, const T2 &x)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != x.size())
		error("illegal shape!");
#endif
	Long i, j, Nr_a = a.nrows(), Nc_a = a.ncols();
	y.resize(a.nrows());
	vecset(y.ptr(), T::value_type(), Nr_a);
	for (i = 0; i < Nr_a; ++i) {
		for (j = 0; j < Nc_a; ++j)
			y[i] += a(i, j) * x[j];
	}
}

template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Matrix<T1> &a, const Vector<T2> &x)
{ mul_vmv(y, a, x); }

template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Cmat<T1> &a, const Vector<T2> &x)
{ mul_vmv(y, a, x); }

// vector-matrix multiplication (row vector assumed)

// parallel version
template <class T, class T1, class T2>
inline void mul_par(Vector<T> &y, const Vector<T1> &x, const Matrix<T2> &a)
{
#ifdef _CHECKBOUNDS_
	if (x.size() != a.nrows()) error("wrong size!");
#endif
	Long j, m{ a.nrows() }, n{ a.ncols() };
	y.resize(n); y = 0.;
	#pragma omp parallel for
	for (j = 0; j < n; ++j) {
		Long k;
		for (k = 0; k < m; ++k)
			y[j] += x[k] * a[k][j];
	}
}

template <class T, class T1, class T2>
inline void mul_vvm(T &y, const T1 &x, const T2 &a)
{
#ifdef _CHECKBOUNDS_
	if (x.size() != a.nrows())
		error("illegal shape!");
#endif
	Long Nr_a = a.nrows(), Nc_a = a.ncols();
	Long i, j, k;
	y.resize(Nc_a);
	vecset(y, T(), Nc_a);
	for (j = 0; j < Nc_a; ++j) {
		for (i = 0; i < Nr_a; ++i)
			y[j] += x[i] * a(i, j);
	}
}

template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Vector<T1> &x, const Matrix<T> &a)
{ mul_vvm(y, x, a); }

template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Vector<T1> &x, const Cmat<T> &a)
{ mul_vvm(y, x, a); }

// matrix-matrix multiplication
// TODO: optimize

template <class T, class T1, class T2>
inline void mul_mmm(T &y, const T1 &a, const T2 &x)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != x.nrows())
		error("illegal shape!");
#endif
	Long Nr_a = a.nrows(), Nc_a = a.ncols(), Nc_x = x.ncols();
	Long i, j, k;
	y.resize(Nr_a, Nc_x);
	vecset(y.ptr(), T::value_type(), Nr_a);
	for (i = 0; i < Nr_a; ++i) {
		for (j = 0; j < Nc_x; ++j) {
			for (k = 0; k < Nc_a; ++k)
				y(i, j) += a(i, k) * x(k, j);
		}
	}
}

template <class T, class T1, class T2>
inline void mul(Matrix<T> &c, const Matrix<T1> &a, const Matrix<T2> &b)
{ mul_mmm(c, a, b); }

template <class T, class T1, class T2>
inline void mul(Cmat<T> &c, const Cmat<T1> &a, const Cmat<T2> &b)
{ mul_mmm(c, a, b); }

// === numerical integration ===

// indefinite integral;
// use cumsum(y)*dx instead
template <class T, class T1>
void cumsum(Vbase<T> &F, const Vbase<T1> &f)
{
	Long i, N{ f.size() };
	F.resize(N); F(0) = f(0);
	for (i = 1; i < N - 1; ++i)
		F(i) = F(i-1) + f(i);
}

// string utilities

template <typename T>
inline std::string num2str(T s)
{
	std::string str = std::to_string(s);
	if (str.find('.') != std::string::npos)
		str.erase(str.find_last_not_of('0') + 1);
	return str;
}

} // namespace slisc
