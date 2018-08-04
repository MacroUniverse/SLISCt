// nr3.h extension
// written by Hongyu Shi
#pragma once
#include "nr3.h"
#include <chrono>

// === constants ===

const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);

// === time utilities ===

extern std::chrono::steady_clock::time_point tic_time_start;
extern std::vector<std::chrono::steady_clock::time_point> tic_time_starts;

inline void tic() { tic_time_start = std::chrono::steady_clock::now(); }

inline Doub toc() {
	std::chrono::steady_clock::time_point tic_time_stop = std::chrono::steady_clock::now();
	std::chrono::duration<double> t = std::chrono::duration_cast<std::chrono::duration<double>>(tic_time_stop - tic_time_start);
	return t.count();
}

inline void tic(Int ind) { tic_time_starts[ind] = std::chrono::steady_clock::now(); }

inline Doub toc(Int ind) {
	std::chrono::steady_clock::time_point tic_time_stop = std::chrono::steady_clock::now();
	std::chrono::duration<double> t = std::chrono::duration_cast<std::chrono::duration<double>>(tic_time_stop - tic_time_starts[ind]);
	return t.count();
}

// === scalar utilities ===

// return 1 if n is odd, return 0 otherwise
inline Int isodd(Int n) { return n & 1; }
inline Long isodd(Long n) { return n & 1; }

// return true if n is power of 2 or 0
inline Bool ispow2(Int n) { return (n&(n-1)) == 0; }
inline Bool ispow2(Long n) { return (n&(n-1)) == 0; }

// operators between Comp and Int
inline Comp operator+(const Comp c, Int_I i) { return c + (Doub)i; }
inline Comp operator+(Int_I i, const Comp c) { return c + (Doub)i; }
inline Comp operator-(Int_I i, const Comp c) { return (Doub)i - c; }
inline Comp operator-(const Comp c, Int_I i) { return c - (Doub)i; }
inline Comp operator*(const Comp c, Int_I i) { return c*(Doub)i; }
inline Comp operator*(Int_I i, const Comp c) { return c*(Doub)i; }
inline Comp operator/(const Comp c, Int_I i) { return c / (Doub)i; }
inline Comp operator/(Int_I i, const Comp c) { return (Doub)i / c; }

// operators between Comp and Long
inline Comp operator+(const Comp c, Long_I i) { return c + (Doub)i; }
inline Comp operator+(Long_I i, const Comp c) { return c + (Doub)i; }
inline Comp operator-(Long_I i, const Comp c) { return (Doub)i - c; }
inline Comp operator-(const Comp c, Long_I i) { return c - (Doub)i; }
inline Comp operator*(const Comp c, Long_I i) { return c*(Doub)i; }
inline Comp operator*(Long_I i, const Comp c) { return c*(Doub)i; }
inline Comp operator/(const Comp c, Long_I i) { return c / (Doub)i; }
inline Comp operator/(Long_I i, const Comp c) { return (Doub)i / c; }

// display vectors and matrices
// don't use template so disp() can be call when debugging
// version 1
void disp(VecUchar_I &v);
void disp(VecInt_I &v);
void disp(VecDoub_I &v);
void disp(VecComp_I &v);
void disp(MatUchar_I &a);
void disp(MatInt_I &a);
void disp(MatDoub_I &a);
void disp(MatComp_I &a);
void disp(Mat3Doub_I &a);
void disp(Mat3Comp_I &a);
// version 2
void disp(VecUchar_I &v, Int_I precision);
void disp(VecInt_I &v, Int_I precision);
void disp(VecDoub_I &v, Int_I precision);
void disp(VecComp_I &v, Int_I precision);
void disp(MatUchar_I &a, Int_I precision);
void disp(MatInt_I &a, Int_I precision);
void disp(MatDoub_I &a, Int_I precision);
void disp(MatComp_I &a, Int_I precision);
void disp(Mat3Doub_I &a, Int_I precision);
void disp(Mat3Comp_I &a, Int_I precision);
// version 3
void disp(VecUchar_I &v, Long_I start, Long_I n);
void disp(VecInt_I &v, Long_I start, Long_I n);
void disp(VecDoub_I &v, Long_I start, Long_I n);
void disp(VecComp_I &v, Long_I start, Long_I n);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
// version 4
void disp(VecUchar_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecInt_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecDoub_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecComp_I &v, Long_I start, Long_I n, Int_I precision);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);

// === get vec/mat properties ===

// total number of elements
template <class T>
inline Long numel(const NRvector<T> &v) { return v.size(); }
template <class T>
inline Long numel(const NRmatrix<T> &a) { return a.nrows()*a.ncols(); }
template <class T>
inline Long numel(const NRmat3d<T> &a) { return a.dim1()*a.dim2()*a.dim3(); }

// pointer to the first element
template <class T>
inline const T* pointer(const NRvector<T> &v) { return &v[0]; }
template <class T>
inline T* pointer(NRvector<T> &v) { return &v[0]; }
template <class T>
inline const T* pointer(const NRmatrix<T> &v) { return v[0]; }
template <class T>
inline T* pointer(NRmatrix<T> &v) { return v[0]; }
template <class T>
inline const T* pointer(const NRmat3d<T> &v) { return v[0][0]; }
template <class T>
inline T* pointer(NRmat3d<T> &v) { return v[0][0]; }

template <class T>
inline T sum(const NRvector<T> &v)
{
	Long i, n{ v.size() };
	T sum = 0;
	for (i = 0; i < n; ++i) {
		sum += v[i];
	}
	return sum;
}

template <class T>
inline T sum(const NRmatrix<T> &a)
{
	Long i, n{ numel(a) };
	T sum = 0;
	const T *pa = a[0];
	for (i = 0; i < n; ++i) {
		sum += pa[i];
	}
	return sum;
}

template <class T>
inline T sum(const NRmat3d<T> &a)
{
	Long i, n{ numel(a) };
	T sum = 0;
	const T *pa = a[0][0];
	for (i = 0; i < n; ++i) {
		sum += pa[i];
	}
	return sum;
}

inline Doub max(VecDoub_I &a)
{
	Long i, N = numel(a);
	Doub val = -1e300;
	const Doub *pa = &a[0];
	for (i = 0; i < N; ++i) {
		if (pa[i] > val) {
			val = pa[i];
		}
	}
	return val;
}

// return max(abs(a(:)))
inline Doub max(VecComp_I &a)
{
	Long i, N = numel(a);
	Doub val = 0., mod;
	const Comp *pa = &a[0];
	for (i = 0; i < N; ++i)
		if ((mod = abs(pa[i])) > val)
			val = mod;
	return val;
}

inline Doub max(Long_O &ind, VecDoub_I &a)
{
	Long i, N = numel(a);
	Doub val = -1e300;
	Doub *pa;
	for (i = 0; i < N; ++i) {
		if (pa[i] > val) {
			val = pa[i];
			ind = i;
		}
	}
	return val;
}

inline Doub max(MatDoub_I &a)
{
	Long i, N = numel(a);
	Doub val = -1e300;
	const Doub *pa = a[0];
	for (i = 0; i < N; ++i) {
		if (pa[i] > val) {
			val = pa[i];
		}
	}
	return val;
}

inline Doub max(Long &ind, MatDoub_I &a)
{
	Long i, N = numel(a);
	Doub val = -1e300;
	Doub *pa;
	for (i = 0; i < N; ++i) {
		if (pa[i] > val) {
			val = pa[i];
			ind = i;
		}
	}
	return val;
}

// norm and norm square of a vector
inline Doub norm(VecDoub_I &v)
{
	Long i, N{ v.size() };
	Doub s{};
	for (i = 0; i < N; ++i) {
		s += v[i] * v[i];
	}
	return sqrt(s);
}

inline Doub norm2(VecDoub_I &v)
{
	Long i, N{ v.size() };
	Doub s{};
	for (i = 0; i < N; ++i) {
		s += v[i] * v[i];
	}
	return s;
}

inline Doub norm(VecComp_I &v)
{
	Long i, N{ v.size() };
	Doub s{}, temp;
	for (i = 0; i < N; ++i) {
		temp = abs(v[i]);
		s += temp*temp;
	}
	return sqrt(s);
}

inline Doub norm2(VecComp_I &v)
{
	Long i, N{ v.size() };
	Doub s{}, temp;
	for (i = 0; i < N; ++i) {
		temp = abs(v[i]);
		s += temp * temp;
	}
	return s;
}

inline Doub norm2(MatComp_I &v)
{
	Long i, N{ numel(v) };
	Doub s{}, temp;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		temp = abs(pv[i]);
		s += temp * temp;
	}
	return s;
}

inline Doub norm2(Mat3Comp_I &v)
{
	Long i, N{ numel(v) };
	Doub s{}, temp;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		temp = abs(pv[i]);
		s += temp * temp;
	}
	return s;
}

// === matrix manipulation ===

// element-wise operators for vectors and matrices

// does not work for integers
template <class T>
inline void linspace(NRvector<T> &v, const T first, const T last, Llong_I n = -1)
{
	Long i, N;
	if (n < 0) N = v.size();
	else v.resize(N = n);	
	T delta = (last - first) / ((T)N - 1);
	for (i = 0; i < N; ++i)
		v[i] = first + delta * (T)i;
}

template <class T>
inline void linspace(NRmatrix<T> &a, const T first, const T last, Llong_I n = -1, Llong_I m = -1)
{
	Long i, N;
	if (n < 0 || m < 0) N = numel(a);
	else {
		a.resize(n, m); N = n*m;
	}
	T delta = (last - first) / ((T)N - 1);
	T *pa = a[0];
	for (i = 0; i < N; ++i)
		pa[i] = first + delta * (T)i;
}

template <class T>
inline void linspace(NRmat3d<T> &a, const T first, const T last, Llong_I n = -1, Llong_I m = -1, Llong_I k = -1)
{
	Long i, N;
	if (n < 0 || m < 0 || k < 0) N = numel(a);
	else {
		a.resize(n, m, k); N = n*m*k;
	}
	T delta = (last - first) / (N - 1);
	T *pa = a[0][0];
	for (i = 0; i < N; ++i)
		pa[i] = first + delta * i;
}

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
inline void flip(NRvector<T> &v)
{
	Long i, n{ v.size() }, ind;
	T temp;
	for (i = 0; i < n / 2; ++i) {
		ind = n - i - 1;
		temp = v[i]; v[i] = v[ind]; v[ind] = temp;
	}
}

template <class T>
inline void flip(NRvector<T> &v, const NRvector<T> &v0)
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
void shift(NRmatrix<T> &a, Llong nshift, Int_I dim = 1)
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
		T *temp = (T*)malloc(sz);
		for (i = 0; i < Nr; ++i) {
			memcpy(&temp[0], a[i], sz);
			memcpy(a[i], a[i] + n, sz_);
			memcpy(a[i] + Nc-n, &temp[0], sz);
		}
		free(temp);
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
		T *temp = (T*)malloc(sz);
		memcpy(temp, a[0], sz);
		memcpy(a[0], a[n], sz_);
		memcpy(a[Nr-n], temp, sz);
		free(temp);
	}
}

// shift the i-th line i times to the left, moving the diagonal to the first column
template <class T>
void diagonals(NRmatrix<T> &a)
{
	Long i, Nr = a.nrows(), Nc = a.ncols();
	T *temp = new T[Nc];
	Long szT = sizeof(T);
	for (i = 1; i < Nr; ++i) {
		memcpy(temp, a[i], i*szT);
		memcpy(a[i], a[i] + i, (Nc-i)*szT);
		memcpy(a[i] + Nc-i, temp, i*szT);
	}
}

template <class T>
void idiagonals(NRmatrix<T> &a)
{
	Long i, Nr = a.nrows(), Nc = a.ncols();
	T *temp = new T[Nc];
	Long szT = sizeof(T);
	for (i = 1; i < Nr; ++i) {
		memcpy(temp, a[i], (Nc-i)*szT);
		memcpy(a[i], a[i] + (Nc-i), i*szT);
		memcpy(a[i] + i, temp, (Nc-i)*szT);
	}
}


// === vectorized math functions ===

template <class T>
void sin(NRvector<T> &y, const NRvector<T> &x)
{
	Long i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = sin(x[i]);
	}
}

template <class T>
void sin(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Long i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = sin(px[i]);
	}
}

template <class T>
void cos(NRvector<T> &y, const NRvector<T> &x)
{
	Long i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = cos(x[i]);
	}
}

template <class T>
void cos(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Long i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = cos(px[i]);
	}
}

template <class T>
void exp(NRvector<T> &y, const NRvector<T> &x)
{
	Long i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = exp(x[i]);
	}
}

template <class T>
void exp(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Long i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = exp(px[i]);
	}
}

template <class T>
void tan(NRvector<T> &y, const NRvector<T> &x)
{
	Long i, N{ x.size() };
	for (i = 0; i < N; ++i) {
		y[i] = tan(x[i]);
	}
}

template <class T>
void tan(NRmatrix<T> &y, const NRmatrix<T> &x)
{
	Long i, N{ numel(x) };
	T *px{ &y[0][0] }, *py{ &x[0][0] };
	for (i = 0; i < N; ++i) {
		py[i] = tan(px[i]);
	}
}

// === matrix arithmatics ===

// operators +=,-=,*=,/= scalar/vec/mat, whenever make sense
template <class T>
inline void operator+=(NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] += v2[i];
}

template <class T>
inline void operator+=(NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] += pv2[i];
}

template <class T>
inline void operator+=(NRmat3d<T> &v1, const NRmat3d<T> &v2)
{
	Long i, N{ v1.dim1()*v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] += pv2[i];
}

template <class T>
inline void operator-=(NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] -= v2[i];
}

template <class T>
inline void operator-=(NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] -= pv2[i];
}

template <class T>
inline void operator-=(NRmat3d<T> &v1, const NRmat3d<T> &v2)
{
	Long i, N{ v1.dim1()*v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] -= pv2[i];
}

template <class T>
inline void operator*=(NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] *= v2[i];
}

template <class T>
inline void operator*=(NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] *= pv2[i];
}

template <class T>
inline void operator*=(NRmat3d<T> &v1, const NRmat3d<T> &v2)
{
	Long i, N{ v1.dim1()*v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] *= pv2[i];
}

template <class T>
inline void operator/=(NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v1[i] /= v2[i];
}

template <class T>
inline void operator/=(NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N{ v1.nrows()*v1.ncols() };
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv1[i] /= pv2[i];
}

template <class T>
inline void operator/=(NRmat3d<T> &v1, const NRmat3d<T> &v2)
{
	Long i, N{ v1.dim1()*v1.dim2()*v1.dim3() };
	auto pv1 = v1[0][0];
	auto pv2 = v2[0][0];
	for (i = 0; i < N; ++i)
		pv1[i] /= pv2[i];
}

template <class T>
inline void operator+=(NRvector<T> &v, Doub_I s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] += s;
}

template <class T>
inline void operator+=(NRmatrix<T> &v, Doub_I s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] += s;
}

template <class T>
inline void operator+=(NRmat3d<T> &v, Doub_I s)
{
	Long i, N{ v.dim1()*v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] += s;
}

template <class T>
inline void operator-=(NRvector<T> &v, Doub_I s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] -= s;
}

template <class T>
inline void operator-=(NRmatrix<T> &v, Doub_I s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] -= s;
}

template <class T>
inline void operator-=(NRmat3d<T> &v, Doub_I s)
{
	Long i, N{ v.dim1()*v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] -= s;
}

template <class T>
inline void operator*=(NRvector<T> &v, Doub_I s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] *= s;
}

template <class T>
inline void operator*=(NRmatrix<T> &v, Doub_I s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] *= s;
}

template <class T>
inline void operator*=(NRmat3d<T> &v, Doub_I s)
{
	Long i, N{ v.dim1()*v.dim2()*v.dim3() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] *= s;
}

template <class T>
inline void operator/=(NRvector<T> &v, Doub_I s)
{
	Long i, N{ v.size() };
	Doub sInv = 1./s;
	for (i = 0; i < N; ++i)
		v[i] *= sInv;
}

template <class T>
inline void operator/=(NRmatrix<T> &v, Doub_I s)
{
	Long i, N{ v.nrows()*v.ncols() };
	Doub sInv = 1./s;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] *= sInv;
}

template <class T>
inline void operator/=(NRmat3d<T> &v, Doub_I s)
{
	Long i, N{ v.dim1()*v.dim2()*v.dim3() };
	Doub sInv = 1./s;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] *= sInv;
}

template <class T>
inline void operator+=(NRvector<T> &v, const Comp s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] += s;
}

template <class T>
inline void operator+=(NRmatrix<T> &v, const Comp s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] += s;
}

template <class T>
inline void operator-=(NRvector<T> &v, const Comp s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] -= s;
}

template <class T>
inline void operator-=(NRmatrix<T> &v, const Comp s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] -= s;
}

template <class T>
inline void operator*=(NRvector<T> &v, const Comp s)
{
	Long i, N{ v.size() };
	for (i = 0; i < N; ++i)
		v[i] *= s;
}

template <class T>
inline void operator*=(NRmatrix<T> &v, const Comp s)
{
	Long i, N{ v.nrows()*v.ncols() };
	auto pv = pointer(v);
	for (i = 0; i < N; ++i)
		pv[i] *= s;
}

template <class T>
inline void operator/=(NRvector<T> &v, const Comp s)
{
	Long i, N{ v.size() };
	Comp sInv = 1./s;
	for (i = 0; i < N; ++i)
		v[i] *= sInv;
}

template <class T>
inline void operator/=(NRmatrix<T> &v, const Comp s)
{
	Long i, N{ v.nrows()*v.ncols() };
	Comp sInv = 1./s;
	auto pv = pointer(v);
	for (i = 0; i < N; ++i) {
		pv[i] *= sInv;
	}
}

template <class T>
inline void plus(NRvector<T> &v, const NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	v.resize(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

template <class T>
inline void plus(NRmatrix<T> &v, const NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N = numel(v1);
	v.resize(v1);
	auto pv = v[0];
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] + pv2[i];
}

template <class T>
inline void plus(NRmat3d<T> &a, const NRmat3d<T> &a1, const NRmat3d<T> &a2)
{
	Long i, N = numel(a1);
	a.resize(a1);
	auto pa = a[0][0];
	auto pa1 = a1[0][0];
	auto pa2 = a2[0][0];
	for (i = 0; i < N; ++i)
		pa[i] = pa1[i] + pa2[i];
}

template <class T>
inline void minus(NRvector<T> &v, const NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	v.resize(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] - v2[i];
}

template <class T>
inline void minus(NRmatrix<T> &v, const NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N = numel(v1);
	v.resize(v1);
	auto pv = v[0];
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] - pv2[i];
}

template <class T>
inline void minus(NRmat3d<T> &a, const NRmat3d<T> &a1, const NRmat3d<T> &a2)
{
	Long i, N = numel(a1);
	a.resize(a1);
	auto pa = a[0][0];
	auto pa1 = a1[0][0];
	auto pa2 = a2[0][0];
	for (i = 0; i < N; ++i)
		pa[i] = pa1[i] - pa2[i];
}

template <class T>
inline void times(NRvector<T> &v, const NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	v.resize(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] * v2[i];
}

template <class T>
inline void times(NRmatrix<T> &v, const NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N = numel(v1);
	v.resize(v1);
	auto pv = v[0];
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] * pv2[i];
}

template <class T>
inline void times(NRmat3d<T> &a, const NRmat3d<T> &a1, const NRmat3d<T> &a2)
{
	Long i, N = numel(a1);
	a.resize(a1);
	auto pa = a[0][0];
	auto pa1 = a1[0][0];
	auto pa2 = a2[0][0];
	for (i = 0; i < N; ++i)
		pa[i] = pa1[i] * pa2[i];
}

template <class T>
inline void divide(NRvector<T> &v, const NRvector<T> &v1, const NRvector<T> &v2)
{
	Long i, N{ v1.size() };
	v.resize(N);
	for (i = 0; i < N; ++i)
		v[i] = v1[i] / v2[i];
}

template <class T>
inline void divide(NRmatrix<T> &v, const NRmatrix<T> &v1, const NRmatrix<T> &v2)
{
	Long i, N = numel(v1);
	v.resize(v1);
	auto pv = v[0];
	auto pv1 = v1[0];
	auto pv2 = v2[0];
	for (i = 0; i < N; ++i)
		pv[i] = pv1[i] / pv2[i];
}

template <class T>
inline void divide(NRmat3d<T> &a, const NRmat3d<T> &a1, const NRmat3d<T> &a2)
{
	Long i, N = numel(a1);
	a.resize(a1);
	auto pa = a[0][0];
	auto pa1 = a1[0][0];
	auto pa2 = a2[0][0];
	for (i = 0; i < N; ++i)
		pa[i] = pa1[i] / pa2[i];
}

inline void real(MatDoub_O &rc, MatComp_I &c)
{
	Long i, m = c.nrows(), n = c.ncols(), N = m*n;
	rc.resize(m, n);
	Doub *prc;
	Comp *pc;
	for (i = 0; i < N; ++i)
		prc[i] = real(pc[i]);
}

inline void imag(MatDoub_O &ic, MatComp_I &c)
{
	Long i, m = c.nrows(), n = c.ncols(), N = m*n;
	ic.resize(m, n);
	Doub *pic;
	Comp *pc;
	for (i = 0; i < N; ++i)
		pic[i] = imag(pc[i]);
}

template <class T>
inline void abs(NRvector<T> &v)
{
	Long i, N = v.size();
	for (i = 0; i < N; ++i)
		v[i] = abs(v[i]);
}

template <class T>
inline void abs(NRmatrix<T> &a)
{
	Long i, N = numel(a);
	const T *p = a[0];
	for (i = 0; i < N; ++i)
		p[i] = abs(p[i]);
}

template <class T>
inline void abs(NRmat3d<T> &a)
{
	Long i, N = numel(a);
	T *p = a[0];
	for (i = 0; i < N; ++i)
		p[i] = abs(p[i]);
}

template <class T, class T1>
inline void abs(NRvector<T1> &out, const NRvector<T> &v)
{
	Long i, N = v.size();
	out.resize(N);
	for (i = 0; i < N; ++i)
		out[i] = abs(v[i]);
}

template <class T, class T1>
inline void abs(NRmatrix<T1> &out, const NRmatrix<T> &a)
{
	Long i, N = numel(a);
	out.resize(a.nrows(), a.ncols());
	T1 *pout = out[0];
	const T *pa = a[0];
	for (i = 0; i < N; ++i)
		pout[i] = abs(pa[i]);
}

template <class T, class T1>
inline void abs(NRmat3d<T1> &out, const NRmat3d<T> &a)
{
	Long i, N = numel(a);
	out.resize(a);
	T1 *pout = out[0];
	T2 *pa = a[0];
	for (i = 0; i < N; ++i)
		pout[i] = abs(pa[i]);
}

inline void complex(VecComp_O &y, VecDoub_I &x)
{
	Long i, N = x.size();
	y.resize(N);
	Doub *py = (Doub*)&y[0];
	for (i = 0; i < N; ++i)
		py[2*i] = x[i];
}

inline void conjugate(VecComp_IO &v)
{
	Long i, N = v.size();
	Doub *p = (Doub *)&v[0];
	N *= 2;
	for (i = 1; i < N; i += 2)
		p[i] = -p[i];
}

// dot products ( conj(v1[i])*v2[i] )
inline Doub operator*(VecDoub_I &v1, VecDoub_I &v2)
{
	Long i, N{ v1.size() };
	Doub s{};
	for (i = 0; i < N; ++i)
		s += v1[i] * v2[i];
	return s;
}

inline Comp operator*(VecComp_I &v1, VecComp_I &v2)
{
	Long i, N{ v1.size() };
	Comp s{};
	for (i = 0; i < N; ++i)
		s += conj(v1[i]) * v2[i];
	return s;
}

inline Comp operator*(VecDoub_I &v1, VecComp_I &v2)
{
	Long i, N{ v1.size() };
	Comp s{};
	for (i = 0; i < N; ++i)
		s += v1[i] * v2[i];
	return s;
}

inline Comp operator*(VecComp_I &v1, VecDoub_I &v2)
{
	Long i, N{ v1.size() };
	Comp s{};
	for (i = 0; i < N; ++i)
		s += conj(v1[i]) * v2[i];
	return s;
}

// outer product ( conj(v1[i})*v2[j] )
inline void outprod(MatComp_O &prod, VecComp_I &v1, VecComp_I &v2)
{
	Long i, j, N1 = v1.size(), N2 = v2.size();
	Comp *pc, v1_i;
	prod.resize(N1, N2);
	for (i = 0; i < N1; ++i) {
		pc = prod[i];
		v1_i = conj(v1[i]);
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}
}

// matrix-vector multiplications (column vector assumed)
template <class T, class T1, class T2>
inline void mul(NRvector<T> &y, const NRmatrix<T1> &a, const NRvector<T2> &x)
{
	Long i, k, m{ a.nrows() }, n{ a.ncols() };
	y.resize(m); y = 0.;
	for (i = 0; i < m; ++i) {
		for (k = 0; k < n; ++k)
			y[i] += a[i][k] * x[k];
	}
}

// vector-matrix multiplication (row vector assumed)
template <class T, class T1, class T2>
inline void mul(NRvector<T> &y, const NRvector<T1> &x, const NRmatrix<T2> &a)
{
	Long j, k, m{ a.nrows() }, n{ a.ncols() };
	y.resize(n); y = 0.;
	for (j = 0; j < n; ++j) {
		for (k = 0; k < m; ++k)
			y[j] += x[k] * a[k][j];
	}
}

// matrix-matrix multiplication
template <class T, class T1, class T2>
inline void mul(NRmatrix<T> &c, const NRmatrix<T1> &a, const NRmatrix<T2> &b)
{
	Long i, j, k, m{ a.nrows() }, n{ b.ncols() }, Nk{ a.ncols() };
	c.resize(m, n); c = 0.;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			for (k = 0; k < Nk; ++k)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
}

// === FT related ===

template <class T>
void fftshift(NRmatrix<T> &a, Int_I dim = 1)
{
	Long m{ a.nrows() }, n{ a.ncols() };
	if (dim == 1) {
		if (isodd(m))
			error("fftshift only supports even rows!")
		else {
			Long halfm = m/2;
			NRmatrix<T> temp(halfm, n);
			Long size = halfm*n*sizeof(T);
			memcpy(temp[0], a[0], size);
			memcpy(a[0], a[halfm], size);
			memcpy(a[halfm], temp[0], size);
		}
	}
	else if (dim == 2) {
		if (isodd(n))
			error("fftshift only supports even columns!")
		else {
			Long i, halfn{ n/2 };
			NRvector<T> temp(halfn);
			Long size{ halfn*sizeof(T) };
			for (i = 0; i < m; ++i) {
				memcpy(&temp[0], a[i], size);
				memcpy(a[i], &a[i][halfn], size);
				memcpy(&a[i][halfn], &temp[0], size);
			}
			
		}
	}
}

// discrete fourier transform from X(x) to Y(k), no fftshift is needed
// each column of X is transformed to each column of Y
// using sum instead of integration, result not normalized
// for each column, Y_j = sum_i ( X_i*exp(-I*k_j*X_i) )
// this is much slower than fft, but for small (xmax-xmin) and (kmax-kmin), could be faster
void dft(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax);

// the inverse of dft, multiplied by 2*pi/(dx*dk).
// this might not be a precise inverse
void idft(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax);
