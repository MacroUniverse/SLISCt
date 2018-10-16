#pragma once
#include "slisc.h"

namespace slisc {

// === scalar utilities ===

// return 1 if n is odd, return 0 otherwise
inline Int isodd(Int_I n) { return n & 1; }
inline Long isodd(Long_I n) { return n & 1; }

// return true if n is power of 2 or 0
inline Bool ispow2(Int_I n) { return (n&(n-1)) == 0; }
inline Bool ispow2(Long_I n) { return (n&(n-1)) == 0; }

// return the positive mod (use "%" when i >= 0)
inline Int mod(Int_I i, Int_I n) { return (i % n + n) % n; }
inline Long_I mod(Long_I i, Long_I n) { return (i % n + n) % n; }

inline Doub sinc(Doub_I x) { return x == 0 ? 1. : std::sin(x)/x; }

// operators between Comp and Int
inline Comp operator+(Comp_I c, Int_I i) { return c + (Doub)i; }
inline Comp operator+(Int_I i, Comp_I c) { return c + (Doub)i; }
inline Comp operator-(Int_I i, Comp_I c) { return (Doub)i - c; }
inline Comp operator-(Comp_I c, Int_I i) { return c - (Doub)i; }
inline Comp operator*(Comp_I c, Int_I i) { return c*(Doub)i; }
inline Comp operator*(Int_I i, Comp_I c) { return c*(Doub)i; }
inline Comp operator/(Comp_I c, Int_I i) { return c / (Doub)i; }
inline Comp operator/(Int_I i, Comp_I c) { return (Doub)i / c; }

// operators between Comp and Long
inline Comp operator+(Comp_I c, Long_I i) { return c + (Doub)i; }
inline Comp operator+(Long_I i, Comp_I c) { return c + (Doub)i; }
inline Comp operator-(Long_I i, Comp_I c) { return (Doub)i - c; }
inline Comp operator-(Comp_I c, Long_I i) { return c - (Doub)i; }
inline Comp operator*(Comp_I c, Long_I i) { return c*(Doub)i; }
inline Comp operator*(Long_I i, Comp_I c) { return c*(Doub)i; }
inline Comp operator/(Comp_I c, Long_I i) { return c / (Doub)i; }
inline Comp operator/(Long_I i, Comp_I c) { return (Doub)i / c; }

// === get vec/mat properties ===

template <class T1, class T2>
inline bool equals_to0(const Vbase<T1> &v1, const Vbase<T2> &v2)
{
	Long i, N{ v1.size() };
	if (N != v2.size()) return false;
	for (i = 0; i < N; ++i)
		if (v1(i) != v2(i))
			return false;
	return true;
}

template <class T1, class T2>
inline bool operator==(const Vector<T1> &v1, const Vector<T2> &v2)
{ return equals_to0(v1, v2); }

template <class T1, class T2>
inline bool operator==(const Matrix<T1> &v1, const Matrix<T2> &v2)
{ return equals_to0(v1, v2); }

template <class T1, class T2>
inline bool operator==(const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ return equals_to0(v1, v2); }

template <class T1, class T2>
inline bool operator!=(const Vector<T1> &v1, const Vector<T2> &v2)
{ return !equals_to0(v1, v2); }

template <class T1, class T2>
inline bool operator!=(const Matrix<T1> &v1, const Matrix<T2> &v2)
{ return !equals_to0(v1, v2); }

template <class T1, class T2>
inline bool operator!=(const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{ return !equals_to0(v1, v2); }

template <class T1, class T2>
inline bool equals_to1(const Vbase<T1> &v, const T2 &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		if (v(i) != s)
			return false;
	return true;
}

template <class T1, class T2>
inline bool operator==(const Vector<T1> &v, const T2 &s)
{ return equals_to1(v, s); }

template <class T1, class T2>
inline bool operator==(const Matrix<T1> &v, const T2 &s)
{ return equals_to1(v, s); }

template <class T1, class T2>
inline bool operator==(const Mat3d<T1> &v, const T2 &s)
{ return equals_to1(v, s); }

template <class T1, class T2>
inline bool operator!=(const Vector<T1> &v, const T2 &s)
{ return !equals_to1(v, s); }

template <class T1, class T2>
inline bool operator!=(const Matrix<T1> &v, const T2 &s)
{ return !equals_to1(v, s); }

template <class T1, class T2>
inline bool operator!=(const Mat3d<T1> &v, const T2 &s)
{ return !equals_to1(v, s); }

template <class T1, class T2>
Bool shape_cmp(const Vector<T1> &v1, const Vector<T2> &v2)
{ return v1.size() == v2.size(); }

template <class T1, class T2>
Bool shape_cmp(const Matrix<T1> &a1, const Matrix<T2> &a2)
{ return (a1.nrows() == a2.nrows()) && (a1.ncols() == a2.ncols()); }

template <class T1, class T2>
Bool shape_cmp(const Mat3d<T1> &a1, const Mat3d<T2> &a2)
{ return (a1.dim1() == a2.dim1()) && (a1.dim2() == a2.dim2()) && (a1.dim3() == a2.dim3()); }

template <class T>
inline T sum(const Vbase<T> &v)
{
	Long i, n{ v.size() };
	T sum = 0;
	for (i = 0; i < n; ++i)
		sum += v(i);
	return sum;
}

template <class T>
inline T max(const Vbase<T> &v)
{
	Long i, N{ v.size() };
	T val{ v(0) };
	for (i = 1; i < N; ++i)
		if (v(i) > val)
			val = v(i);
	return val;
}

// for Comp, return max(abs(a(:))
inline Doub max(const Vbase<Comp> &v)
{
	Long i, N{ v.size() };
	Doub val{ abs(v(0)) };
	for (i = 1; i < N; ++i)
		if (abs(v(i)) > val)
			val = abs(v(i));
	return val;
}

template <class T>
inline T max(Long_O &ind, const Vbase<T> &v)
{
	Long i, N{ v.size() };
	T val{ v(0) };
	for (i = 1; i < N; ++i)
		if (v(i) > val) {
			val = v(i); ind = i;
		}
	return val;
}

inline Doub max(Long_O &ind, const Vbase<Comp> &v)
{
	Long i, N{ v.size() };
	Doub val{ abs(v(0)) };
	for (i = 1; i < N; ++i)
		if (abs(v(i)) > val) {
			val = abs(v(i)); ind = i;
		}
	return val;
}

// sum(v(:).^2) for real numbers
template <class T>
inline T norm2(Vbase<T> &v)
{
	Long i, N{ v.size() };
	T s2{};
	for (i = 0; i < N; ++i)
		s2 += SQR(v(i));
	return s2;
}

template <class T>
inline T norm(Vbase<T> &v)
{ return std::sqrt(norm2(v)); }

//sum(abs(v(:)). ^ 2) for complex numbers
inline Doub norm2(Vbase<Comp> &v)
{
	Long i, N{ v.size() };
	Doub s2{};
	for (i = 0; i < N; ++i)
		s2 += SQR(abs(v(i)));
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
		v(i) = first + delta * T(i);
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
inline void her(MatComp_O &h, MatComp_I &a)
{
	Long i, j, m = a.nrows(), n = a.ncols();
	h.resize(n, m);
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			h[j][i] = conj(a[i][j]);
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
inline void sin0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = sin(v1(i));
}

template <class T, class T1>
void sin(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); sin0(v, v1); }

template <class T, class T1>
void sin(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); sin0(v, v1); }

template <class T, class T1>
void sin(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); sin0(v, v1); }

template <class T, class T1>
inline void cos0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) = cos(v1(i));
}

template <class T, class T1>
void cos(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); cos0(v, v1); }

template <class T, class T1>
void cos(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); cos0(v, v1); }

template <class T, class T1>
void cos(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); cos0(v, v1); }

template <class T, class T1>
inline void exp0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) = exp(v1(i));
}

template <class T, class T1>
void exp(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); exp0(v, v1); }

template <class T, class T1>
void exp(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); exp0(v, v1); }

template <class T, class T1>
void exp(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); exp0(v, v1); }

template <class T, class T1>
inline void tan0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) = tan(v1(i));
}

template <class T, class T1>
void tan(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); tan0(v, v1); }

template <class T, class T1>
void tan(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); tan0(v, v1); }

template <class T, class T1>
void tan(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); tan0(v, v1); }

// === matrix arithmatics ===

// v += v
template <class T, class T1>
inline void plus_equals0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) += v1(i);
}

template <class T, class T1>
inline void operator+=(Vector<T> &v, const Vector<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	plus_equals0(v, v1);
}

template <class T, class T1>
inline void operator+=(Matrix<T> &v, const Matrix<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	plus_equals0(v, v1);
}

template <class T, class T1>
inline void operator+=(Mat3d<T> &v, const Mat3d<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	plus_equals0(v, v1);
}

// v -= v
template <class T, class T1>
inline void minus_equals0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) -= v1(i);
}

template <class T, class T1>
inline void operator-=(Vector<T> &v, const Vector<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	minus_equals0(v, v1);
}

template <class T, class T1>
inline void operator-=(Matrix<T> &v, const Matrix<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	minus_equals0(v, v1);
}

template <class T, class T1>
inline void operator-=(Mat3d<T> &v, const Mat3d<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	minus_equals0(v, v1);
}

// v *= v
template <class T, class T1>
inline void times_equals0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) *= v1(i);
}

template <class T, class T1>
inline void operator*=(Vector<T> &v, const Vector<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	times_equals0(v, v1);
}

template <class T, class T1>
inline void operator*=(Matrix<T> &v, const Matrix<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	times_equals0(v, v1);
}

template <class T, class T1>
inline void operator*=(Mat3d<T> &v, const Mat3d<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	times_equals0(v, v1);
}

// v /= v
template <class T, class T1>
inline void divide_equals0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) /= v1(i);
}

template <class T, class T1>
inline void operator/=(Vector<T> &v, const Vector<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	divide_equals0(v, v1);
}

template <class T, class T1>
inline void operator/=(Matrix<T> &v, const Matrix<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	divide_equals0(v, v1);
}

template <class T, class T1>
inline void operator/=(Mat3d<T> &v, const Mat3d<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	divide_equals0(v, v1);
}

template <class T, class T1>
inline void plus_equals1(Vbase<T> &v, const T1 &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) += s;
}

// v += s
template <class T, class T1>
inline void operator+=(Vector<T> &v, const T1 &s)
{ plus_equals1(v, s); }

template <class T, class T1>
inline void operator+=(Matrix<T> &v, const T1 &s)
{ plus_equals1(v, s); }

template <class T, class T1>
inline void operator+=(Mat3d<T> &v, const T1 &s)
{ plus_equals1(v, s); }

// v -= s
template <class T, class T1>
inline void minus_equals1(Vbase<T> &v, const T1 &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) -= s;
}

template <class T, class T1>
inline void operator-=(Vector<T> &v, const T1 &s)
{ minus_equals1(v, s); }

template <class T, class T1>
inline void operator-=(Matrix<T> &v, const T1 &s)
{ minus_equals1(v, s); }

template <class T, class T1>
inline void operator-=(Mat3d<T> &v, const T1 &s)
{ minus_equals1(v, s); }

// v *= s
template <class T, class T1>
inline void times_equals1(Vbase<T> &v, const T1 &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) *= s;
}

template <class T, class T1>
inline void operator*=(Vector<T> &v, const T1 &s)
{ times_equals1(v, s); }

template <class T, class T1>
inline void operator*=(Matrix<T> &v, const T1 &s)
{ times_equals1(v, s); }

template <class T, class T1>
inline void operator*=(Mat3d<T> &v, const T1 &s)
{ times_equals1(v, s); }

// v /= s (only works for floating point types)
template <class T, class T1>
inline void operator/=(Vector<T> &v, const T1 &s)
{ v *= 1./s; }

template <class T, class T1>
inline void operator/=(Matrix<T> &v, const T1 &s)
{ v *= 1./s; }

template <class T, class T1>
inline void operator/=(Mat3d<T> &v, const T1 &s)
{ v *= 1./s; }

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
inline void rem0(Vbase<T> &v, const Vbase<T> &v1, const T &s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) % s;
}

template <class T>
inline void rem(Vector<T> &v, const Vector<T> &v1, const T &s)
{ v.resize(v1); rem0(v, v1, s); }

template <class T>
inline void rem(Matrix<T> &v, const Matrix<T> &v1, const T &s)
{ v.resize(v1); rem0(v, v1, s); }

template <class T>
inline void rem(Mat3d<T> &v, const Mat3d<T> &v1, const T &s)
{ v.resize(v1); rem0(v, v1, s); }

// TODO : rem(v, s, v1)
// TODO : rem(v, v1, v2)

// mod(v, v, s)
template <class T>
inline void mod0(Vbase<T> &v, const Vbase<T> &v1, const T &s)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = mod(v1(i), s);
}

template <class T>
inline void mod(Vector<T> &v, const Vector<T> &v1, const T &s)
{ v.resize(v1); mod0(v, v1, s); }

template <class T>
inline void mod(Matrix<T> &v, const Matrix<T> &v1, const T &s)
{ v.resize(v1); mod0(v, v1, s); }

template <class T>
inline void mod(Mat3d<T> &v, const Mat3d<T> &v1, const T &s)
{ v.resize(v1); mod0(v, v1, s); }

// TODO : mod(v, s, v1)
// TODO : mod(v, v1, v2)

// plus(v, v, s)
template <class T, class T1, class T2>
inline void plus0(Vbase<T> &v, const Vbase<T1> &v1, const T2 &s)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) + s;
}

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ v.resize(v1); plus0(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ v.resize(v1); plus0(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ v.resize(v1); plus0(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ plus(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ plus(v, v1, s); }

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ plus(v, v1, s); }

// plus(v, v, v)
template <class T, class T1, class T2>
inline void plus1(Vbase<T> &v, const Vbase<T1> &v1, const Vbase<T2> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) + v2(i);
}

template <class T, class T1, class T2>
inline void plus(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); plus1(v, v1, v2);
}

template <class T, class T1, class T2>
inline void plus(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); plus1(v, v1, v2);
}

template <class T, class T1, class T2>
inline void plus(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); plus1(v, v1, v2);
}

// minus(v)
template <class T>
inline void minus(Vbase<T> &v)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v(i) *= -1;
}

// minus(v, v)
template <class T, class T1>
inline void minus1(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = -v1(i);
}

template <class T, class T1>
inline void minus(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); minus1(v, v1); }

template <class T, class T1>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); minus1(v, v1); }

template <class T, class T1>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); minus1(v, v1); }

// minus(v, s, v)
template <class T, class T1, class T2>
inline void minus2(Vbase<T> &v, const T1 &s, const Vbase<T2> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = s - v1(i);
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ v.resize(v1); minus2(v, s, v1); }

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ v.resize(v1); minus2(v, s, v1); }

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ v.resize(v1); minus2(v, s, v1); }

// minus(v, v, s)
template <class T, class T1, class T2>
inline void minus3(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{
	Long i, N{ v1.size() };
	v.resize(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ v.resize(v1); minus3(v, v1, s); }

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ v.resize(v1); minus3(v, v1, s); }

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ v.resize(v1); minus3(v, v1, s); }

// minus(v, v, v)
template <class T, class T1, class T2>
inline void minus4(Vbase<T> &v, const Vbase<T1> &v1, const Vbase<T2> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) - v2(i);
}

template <class T, class T1, class T2>
inline void minus(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); minus4(v, v1, v2);
}

template <class T, class T1, class T2>
inline void minus(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); minus4(v, v1, v2);
}

template <class T, class T1, class T2>
inline void minus(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); minus4(v, v1, v2);
}

// times(v, v, s)
template <class T, class T1, class T2>
inline void times0(Vbase<T> &v, const Vbase<T1> &v1, const T2 &s)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) * s;
}

template <class T, class T1, class T2>
inline void times(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ v.resize(v1); times0(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ v.resize(v1); times0(v, v1, s); }

template <class T, class T1, class T2>
inline void times(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ v.resize(v1); times0(v, v1, s); }

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
inline void times1(Vbase<T> &v, const Vbase<T1> &v1, const Vbase<T2> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) * v2(i);
}

template <class T, class T1, class T2>
inline void times(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); times1(v, v1, v2);
}

template <class T, class T1, class T2>
inline void times(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); times1(v, v1, v2);
}

template <class T, class T1, class T2>
inline void times(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); times1(v, v1, v2);
}

// times(v, v, s)
template <class T, class T1, class T2>
inline void divide0(Vbase<T> &v, const Vbase<T1> &v1, const T2 &s)
{
	Long i, N{ v1.size() };
	T2 sInv{ 1/s };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) * sInv;
}

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const Vector<T1> &v1, const T2 &s)
{ v.resize(v1); divide0(v, v1, s); }

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const Matrix<T1> &v1, const T2 &s)
{ v.resize(v1); divide0(v, v1, s); }

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const Mat3d<T1> &v1, const T2 &s)
{ v.resize(v1); divide0(v, v1, s); }

// divide(v, s, v)
template <class T, class T1, class T2>
inline void divide1(Vbase<T> &v, const T1 &s, const Vbase<T2> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = s / v1(i);
}

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const T1 &s, const Vector<T2> &v1)
{ v.resize(v1); divide1(v, s, v1); }

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const T1 &s, const Matrix<T2> &v1)
{ v.resize(v1); divide1(v, s, v1); }

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const T1 &s, const Mat3d<T2> &v1)
{ v.resize(v1); divide1(v, s, v1); }

// divide(v, v, v)
template <class T, class T1, class T2>
inline void divide3(Vbase<T> &v, const Vbase<T1> &v1, const Vbase<T2> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i) / v2(i);
}

template <class T, class T1, class T2>
inline void divide(Vector<T> &v, const Vector<T1> &v1, const Vector<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); divide3(v, v1, v2);
}

template <class T, class T1, class T2>
inline void divide(Matrix<T> &v, const Matrix<T1> &v1, const Matrix<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); divide3(v, v1, v2);
}

template <class T, class T1, class T2>
inline void divide(Mat3d<T> &v, const Mat3d<T1> &v1, const Mat3d<T2> &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); divide3(v, v1, v2);
}

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
inline void real0(Vbase<T> &v, const Vbase<Comp> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = real(v1(i));
}

template <class T>
inline void real(Vector<T> &v, const Vector<Comp> &v1)
{ v.resize(v1); real0(v, v1); }

template <class T>
inline void real(Matrix<T> &v, const Matrix<Comp> &v1)
{ v.resize(v1); real0(v, v1); }

template <class T>
inline void real(Mat3d<T> &v, const Mat3d<Comp> &v1)
{ v.resize(v1); real0(v, v1); }

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
inline void imag0(Vbase<T> &v, const Vbase<Comp> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = imag(v1(i));
}

template <class T>
inline void imag(Vector<T> &v, const Vector<Comp> &v1)
{ v.resize(v1); imag0(v, v1); }

template <class T>
inline void imag(Matrix<T> &v, const Matrix<Comp> &v1)
{ v.resize(v1); imag0(v, v1); }

template <class T>
inline void imag(Mat3d<T> &v, const Mat3d<Comp> &v1)
{ v.resize(v1); imag0(v, v1); }

// abs(v)
template <class T>
inline void abs(Vbase<T> &v)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; i += 2)
		v(i) = std::abs(v(i));
}

// abs(v, v)
template <class T, class T1>
inline void abs0(Vbase<T> &v, const Vbase<T1> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = abs(v1(i));
}

template <class T, class T1>
inline void abs(Vector<T> &v, const Vector<T1> &v1)
{ v.resize(v1); abs0(v, v1); }

template <class T, class T1>
inline void abs(Matrix<T> &v, const Matrix<T1> &v1)
{ v.resize(v1); abs0(v, v1); }

template <class T, class T1>
inline void abs(Mat3d<T> &v, const Mat3d<T1> &v1)
{ v.resize(v1); abs0(v, v1); }

// doubl2comp(v, v)
inline void doub2comp0(Vbase<Comp> &v, const Vbase<Doub> &v1)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v(i) = v1(i);
}

inline void doub2comp(Vector<Comp> &v, const Vector<Doub> &v1)
{ v.resize(v1); doub2comp0(v, v1); }

inline void doub2comp(Matrix<Comp> &v, const Matrix<Doub> &v1)
{ v.resize(v1); doub2comp0(v, v1); }

inline void doub2comp(Mat3d<Comp> &v, const Mat3d<Doub> &v1)
{ v.resize(v1); doub2comp0(v, v1); }

// conj(v)
inline void conj(Vbase<Comp> &v)
{
	Long i, N{ 2 * v.size() };
	Doub *p = (Doub *)v.ptr();
	for (i = 1; i < N; i += 2)
		p[i] = -p[i];
}

// dot products ( sum conj(v1[i])*v2[i] )
// s = dot(v, v)
template <class T, class T1, class T2>
inline T dot0(const Vector<T1> &v1, const Vector<T2> &v2)
{
	Long i, N{ v1.size() };
	T s{};
	for (i = 0; i < N; ++i)
		s += v1[i] * v2[i];
	return s;
}

template <class T, class T1, class T2>
inline T dot1(const Vector<T1> &v1, const Vector<T2> &v2)
{
	Long i, N{ v1.size() };
	T s{};
	for (i = 0; i < N; ++i)
		s += conj(v1[i]) * v2[i];
	return s;
}

inline Doub operator*(VecDoub_I &v1, VecDoub_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	return dot0<Doub>(v1, v2);
}

inline Comp operator*(VecComp_I &v1, VecComp_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	return dot1<Comp>(v1, v2);
}

inline Comp operator*(VecDoub_I &v1, VecComp_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	return dot0<Comp>(v1, v2);
}

inline Comp operator*(VecComp_I &v1, VecDoub_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	return dot1<Comp>(v1, v2);
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
inline void outprod(Matrix<T> &v, VecComp_I &v1, const Vector<T2> &v2)
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
inline void outprod_par(Matrix<T> &v, VecComp_I &v1, const Vector<T2> &v2)
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

// matrix-vector multiplications (column vector assumed)
template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Matrix<T1> &a, const Vector<T2> &x)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != x.size()) error("wrong shape!");
#endif
	Long i, k, m{ a.nrows() }, n{ a.ncols() };
	y.resize(m); y = 0.;
	for (i = 0; i < m; ++i) {
		for (k = 0; k < n; ++k)
			y[i] += a[i][k] * x[k];
	}
}

// vector-matrix multiplication (row vector assumed)
template <class T, class T1, class T2>
inline void mul(Vector<T> &y, const Vector<T1> &x, const Matrix<T2> &a)
{
#ifdef _CHECKBOUNDS_
	if (x.size() != a.nrows()) error("wrong size!");
#endif
	Long j, k, m{ a.nrows() }, n{ a.ncols() };
	y.resize(n); y = 0.;
	for (j = 0; j < n; ++j) {
		for (k = 0; k < m; ++k)
			y[j] += x[k] * a[k][j];
	}
}

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

// matrix-matrix multiplication
// TODO: optimize
template <class T, class T1, class T2>
inline void mul(Matrix<T> &c, const Matrix<T1> &a, const Matrix<T2> &b)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != b.nrows()) error("wrong size!");
#endif
	Long i, j, k, m{ a.nrows() }, n{ b.ncols() }, Nk{ a.ncols() };
	c.resize(m, n); c = 0.;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			for (k = 0; k < Nk; ++k)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
}

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
