#pragma once
#include "global.h"
#include "meta.h"
#include "vector.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "fixsize.h"
#include "ptr_arith.h"

namespace slisc {

// === get vec/mat properties ===

// check if vec/mat sizes are the same
template <class T1, class T2, SLS_IF(is_contain<T1>() && is_contain<T2>())>
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

template <class T1, class T2, SLS_IF(is_dense<T1>() && is_same_contain<T1,T2>())>
Bool operator==(const T1 &v1, const T2 &v2)
{
	return shape_cmp(v1, v2) && equals_to_vv(v1.ptr(), v2.ptr(), v2.size());
}

template <class T1, class T2, SLS_IF(
	is_dense_mat<T1>() && is_dense_mat<T2>() && !is_same_contain<T1, T2>()
)>
Bool operator==(const T1 &v1, const T2 &v2)
{
	if (!shape_cmp(v1, v2)) return false;
	for (Long i = 0; i < v1.nrows(); ++i)
		for (Long j = 0; j < v1.ncols(); ++j)
			if (v1(i, j) != v2(i, j))
				return false;
	return true;
}

template <class Tv, class Ts, SLS_IF(is_dense<Tv>() && !is_contain<Ts>())>
Bool operator==(const Tv &v, const Ts &s)
{
	return equals_to_vs(v.ptr(), s, v.size());
}

template <class Tv, class Ts, SLS_IF(is_dense<Tv>() && !is_contain<Ts>())>
Bool operator==(const Ts &s, const Tv &v)
{ return v == s; }

// operator!= for slisc containers
template <class T1, class T2, SLS_IF(is_dense<T1>() || is_dense<T2>())>
Bool operator!=(const T1 &v1, const T2 &v2)
{
	return !(v1 == v2);
}

// sum of container elements

template <class T, SLS_IF(is_dense<T>())>
inline const auto sum(const T &v)
{
	return sum_v(v.ptr(), v.size());
}

template <class T, SLS_IF(is_dense<T>())>
inline const auto max(const T &v)
{
	return max_v(v.ptr(), v.size());
}

// return max(abs(a(:))
template <class T, SLS_IF(is_dense<T>())>
inline const auto max_abs(const T &v)
{ return max_abs_v(v.ptr(), v.size()); }

template <class T, SLS_IF(is_dense<T>())>
inline const auto max(Long_O ind, const T &v)
{
	Long i, N{ v.size() };
	auto val = v[0];
	for (i = 1; i < N; ++i) {
		if (val < v[i]) {
			val = v[i]; ind = i;
		}
	}
	return val;
}

// s = norm2(v)  (|v1|^2 + |v2|^2 + ...)
template <class T, SLS_IF(is_dense<T>())>
const auto norm2(T &v)
{
	Long i, N{ v.size() };
	auto s2 = ABS2(v[0]);
	for (i = 1; i < N; ++i)
		s2 += ABS2(v[i]);
	return s2;
}

template <class T>
const auto norm(T &v, SLS_IF(is_dense<T>()))
{ return sqrt(norm2(v)); }

// === matrix manipulation ===

// does not work for integers

template <class T, class T1, class T2, SLS_IF(is_scalar<T>())>
inline void linspace(Vector<T> &v, const T1 &first, const T2 &last, Llong_I N = -1)
{
	if (N >= 0) v.resize(N);
	linspace_vss(v.ptr(), (T)first, (T)last, v.size());
}

template <class T, class T1, class T2, SLS_IF(is_dense_mat<T>())>
inline void linspace(T &v, const T1 &first, const T2 &last,
	Llong_I Nr = -1, Long_I Nc = -1)
{
	if (Nr >= 0 && Nc >= 0) v.resize(Nr, Nc);
	linspace_vss(v.ptr(), (contain_type<T>)first, (contain_type<T>)last, v.size());
}

template <class T, class T1, class T2, SLS_IF(is_Mat3d<T>())>
inline void linspace(T &v, const T1 &first, const T2 &last,
	Llong_I N1 = -1, Long_I N2 = -1, Long_I N3 = -1)
{
	if (N1 >= 0 && N2 >= 0 && N3 >= 0) v.resize(N1, N2, N3);
	linspace_vss(v.ptr(), (contain_type<T>)first, (contain_type<T>)last, v.size());
}

// === vectorized math functions ===

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T,T1>())>
void sqrt(T &v, const T1 &v1)
{ v.resize(v1); sqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void invSqrt(T &v, const T1 &v1)
{ v.resize(v1); invSqrt_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void sin(T &v, const T1 &v1)
{ v.resize(v1); sin_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void cos(T &v, const T1 &v1)
{ v.resize(v1); cos_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void exp(T &v, const T1 &v1)
{ v.resize(v1); exp_vv(v.ptr(), v1.ptr(), v1.size()); }

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void tan(T &v, const T1 &v1)
{ v.resize(v1); tan_vv(v.ptr(), v1.ptr(), v1.size()); }

// === vector/matrix arithmetics ===

// flip(v)
template <class T, SLS_IF(is_Vector<T>())>
inline void flip(T &v)
{
	flip_v(v.ptr(), v.size());
}

// v = flip(v)
template <class T, class T1, SLS_IF(
	is_Vector<T>() && is_Vector<T1>()
)>
inline void flip(T &v, const T1 &v1)
{
	v.resize(v1);
	flip(v.ptr(), v1.ptr(), v1.size());
}

// matrix transpose

// trans(v)
template <class T, SLS_IF(is_dense_mat<T>())>
inline void trans(T &v)
{
#ifdef SLS_CHECK_SHAPE
	if (v.nrows() != v.ncols()) error("illegal shape!");
#endif
	for (Long i = 0; i < v.nrows(); ++i)
		for (Long j = 0; j < i; ++j)
			swap(v(i, j), v(j, i));
}

// v = trans(v)
template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>()
)>
inline void trans(T &v, const T1 &v1)
{
	v.resize(v1.ncols(), v1.nrows());
	for (Long i = 0; i < v.nrows(); ++i)
		for (Long j = 0; j < v.ncols(); ++j)
			v(i, j) = v1(j, i);
}

// hermitian conjugate

// her(v)
template <class T, SLS_IF(is_comp_dense<T>() && is_dense_mat<T>())>
inline void her(T &v)
{
#ifdef SLS_CHECK_SHAPE
	if (v.nrows() != v.ncols()) error("illegal shape!");
#endif
	// diagonal
	for (Long i = 0; i < v.nrows(); ++i)
		v(i, i) = conj(v(i, i));
	// off-diagonal
	for (Long i = 0; i < v.nrows(); ++i) {
		for (Long j = 0; j < i; ++j) {
			contain_type<T> temp = v(i, j);
			v(i, j) = conj(v(j, i));
			v(j, i) = conj(temp);
		}
	}
}

// v = her(v)
template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>() &&
	is_comp_dense<T1>() &&
	type_num<contain_type<T>>() >= type_num<contain_type<T1>>()
)>
inline void her(T &v, const T1 &v1)
{
	v.resize(v1.ncols(), v1.nrows());
	for (Long i = 0; i < v.nrows(); ++i) {
		for (Long j = 0; j < v.ncols(); ++j)
			v(i, j) = conj(v1(j, i));
	}
}

// v += v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator+=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	plus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v -= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator-=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	minus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v *= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator*=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	times_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v /= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator/=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	divide_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v += s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator+=(T &v, const Ts &s)
{
	plus_equals_vs(v.ptr(), s, v.size());
}

// v -= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator-=(T &v, const Ts &s)
{
	minus_equals_vs(v.ptr(), s, v.size());
}

// v *= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator*=(T &v, const Ts &s)
{
	times_equals_vs(v.ptr(), s, v.size());
}

// v /= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator/=(T &v, const Ts &s)
{
	divide_equals_vs(v.ptr(), s, v.size());
}

// v %= s
template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void rem(T &v, const Ts &s)
{
	rem_vs(v.ptr(), s, v.size());
}

// v = v % s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T, T1>() && is_scalar<Ts>()
)>
inline void rem(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1);
	rem_vvs(v.ptr(), v1.ptr(), s, v.size());
}

// mod(v, s)
template <class T, class T1, SLS_IF(
	is_dense<T>() && is_scalar<T1>()
)>
inline void mod(T &v, const T1 &s)
{
	mod_vs(v.ptr(), s, v.size());
}

// v = mod(v, s)
template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T, T1>()
)>
inline void mod(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1);
	mod_vvs(v.ptr(), v1.ptr(), s, v.size());
}

// TODO : mod(v, s, v1)
// TODO : mod(v, v1, v2)

// v = v + s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Plus(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1); plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Plus(T &v, const Ts &s, const T1 &v1)
{
	v.resize(v1); plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v + v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Plus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); plus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// -v inplace
template <class T, SLS_IF(is_dense<T>())>
inline void Minus(T &v)
{
	minus_v(v.ptr(), v.size());
}

// v = -v

template <class T, class T1, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>()
)>
inline void Minus(T &v, const T1 &v1)
{
	v.resize(v1); minus_vv(v.ptr(), v1.ptr(), v1.size());
}

// v = s - v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Minus(T &v, const Ts &s, const T1 &v1)
{
	v.resize(v1); minus_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v - s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Minus(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1); minus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v - v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Minus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); minus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

// v = v * s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Times(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1); times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s * v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Times(T &v, const Ts &s, const T1 &v1)
{
	v.resize(v1); times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v * v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Times(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); times_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// v = v / s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Divide(T &v, const T1 &v1, const Ts &s)
{
	v.resize(v1); divide_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s / v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Divide(T &v, const Ts &s, const T1 &v1)
{
	v.resize(v1); divide_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v / v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Divide(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	v.resize(v1); divide_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// real(v)

template <class T, SLS_IF(is_comp_dense<T>())>
inline void real(T &v)
{
	real_v(v.ptr(), v.size());
}

// v = real(v)

template <class T, class T1,
	SLS_IF(is_same_contain<T, T1>(), is_comp_dense<T1>())>
inline void real(T &v, const T1 &v1)
{
	v.resize(v1); real_vv(v.ptr(), v1.ptr(), v1.size());
}

// imag(v)
template <class T, SLS_IF(is_comp_dense<T>())>
inline void imag(T &v)
{
	imag_v(v.ptr(), v.size());
}

// v = imag(v)
template <class T, class T1,
	SLS_IF(is_same_contain<T, T1>() && is_comp_dense<T1>())>
	inline void imag(T &v, const T1 &v1)
{
	v.resize(v1); imag_vv(v.ptr(), v1.ptr(), v1.size());
}

// abs(v)
template <class T, SLS_IF(is_dense<T>())>
inline void abs(T &v)
{
	abs_v(v.ptr(), v.size());
}

// v = abs(v)
template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void abs(T &v, const T1 &v1)
{
	v.resize(v1); abs_vv(v.ptr(), v1.ptr(), v1.size());
}

// to_comp() deleted, use operator= instead

// conj(v)
template <class T, SLS_IF(is_comp_dense<T>())>
inline void conj(T &v)
{
	conj_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(
	is_comp_dense<T>() && is_comp_dense<T1>() &&
	is_same_contain<T,T1>()
)>
inline void conj(T &v, const T1 &v1)
{
	v.resize(v1);
	conj_vv(v.ptr(), v1.ptr(), v1.size());
}

// dot products ( sum conj(v1[i])*v2[i] )

// s = dot(v, v)
template <class T1, class T2, SLS_IF(is_Vector<T1>() && is_Vector<T2>())>
inline auto dot(const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2)) error("wrong shape!");
#endif
	return dot_vv(v1.ptr(), v2.ptr(), v2.size());
}

// outer product ( v(i,j) = conj(v1[i])*v2[j] )
// outprod(v, v, v)
template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_Vector<T1>() && is_Vector<T2>()
)>
inline void outprod(T &v, const T1 &v1, const T2 &v2)
{
	error("TODO");
	/*Long i, j, N1{ v1.size() }, N2{ v2.size() };
	Comp *pc, v1_i;
	v.resize(N1, N2);
	for (i = 0; i < N1; ++i) {
		pc = v[i];
		v1_i = v1[i];
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}*/
}

// parallel version
template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_Vector<T1>() && is_Vector<T2>()
)>
inline void outprod_par(T &v, const T1 &v1, const T2 &v2)
{
	error("TODO");
	/*Long i, N1{ v1.size() }, N2{ v2.size() };
	v.resize(N1, N2);
	#pragma omp parallel for
	for (i = 0; i < N1; ++i) {
		Long j;
		Comp *pc, v1_i;
		pc = v[i];
		v1_i = v1[i];
		for (j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}*/
}

template <class T, class T1, class T2, SLS_IF(
	is_Vector<T>() && is_dense_mat<T1>() && is_Vector<T2>()
)>
inline void mul(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != x.size())
		error("illegal shape!");
#endif
	Long i, j, Nr_a = a.nrows(), Nc_a = a.ncols();
	y.resize(a.nrows());
	vecset(y.ptr(), contain_type<T>(), Nr_a);
	for (i = 0; i < Nr_a; ++i) {
		for (j = 0; j < Nc_a; ++j)
			y[i] += a(i, j) * x[j];
	}
}

// vector-matrix multiplication (row vector assumed)

// parallel version
template <class T, class T1, class T2, SLS_IF(
	is_Vector<T>() && is_Vector<T1>() && is_dense_mat<T2>()
)>
inline void mul_par(T &y, const T1 &x, const T2 &a)
{
	error("TODO");
//#ifdef SLS_CHECK_BOUNDS
//	if (x.size() != a.nrows()) error("wrong size!");
//#endif
//	Long j, m{ a.nrows() }, n{ a.ncols() };
//	y.resize(n); y = 0.;
//	#pragma omp parallel for
//	for (j = 0; j < n; ++j) {
//		Long k;
//		for (k = 0; k < m; ++k)
//			y[j] += x[k] * a[k][j];
//	}
}

template <class T, class T1, class T2, SLS_IF(
	is_Vector<T>() && is_Vector<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &x, const T2 &a)
{
#ifdef SLS_CHECK_SHAPE
	if (x.size() != a.nrows())
		error("illegal shape!");
#endif
	Long Nr_a = a.nrows(), Nc_a = a.ncols();
	Long i, j, k;
	y.resize(Nc_a);
	vecset(y.ptr(), contain_type<T>(), Nc_a);
	for (j = 0; j < Nc_a; ++j) {
		for (i = 0; i < Nr_a; ++i)
			y[j] += x[i] * a(i, j);
	}
}

// matrix-matrix multiplication

template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != x.nrows())
		error("illegal shape!");
#endif
	Long Nr_a = a.nrows(), Nc_a = a.ncols(), Nc_x = x.ncols();
	Long i, j, k;
	y.resize(Nr_a, Nc_x);
	vecset(y.ptr(), contain_type<T>(), Nr_a*Nc_x);
	for (i = 0; i < Nr_a; ++i) {
		for (j = 0; j < Nc_x; ++j) {
			for (k = 0; k < Nc_a; ++k)
				y(i, j) += a(i, k) * x(k, j);
		}
	}
}

// === numerical integration ===

// indefinite integral;
// use cumsum(y)*dx instead
template <class T, class T1, SLS_IF(is_Vector<T>() && is_Vector<T1>())>
void cumsum(T &v, const T1 &v1)
{
	v.resize(v1); cumsum_vv(v.ptr(), v1.ptr(), v1.size());
}

} // namespace slisc
