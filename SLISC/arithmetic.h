#pragma once
#include "global.h"
#include "meta.h"
#include "vector.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "cmat3d.h"
#include "fixsize.h"
#include "ptr_arith.h"

namespace slisc {

// === get vec/mat properties ===

// check if vec/mat sizes are the same
template <class T1, class T2, SLS_IF(is_contain<T1>() && is_contain<T2>())>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
	if constexpr (ndims<T1>() == 1 && ndims<T2>() == 1) {
		return v1.size() == v2.size();
	}
	else if constexpr (ndims<T1>() == 2 && ndims<T2>() == 2) {
		return v1.n1() == v2.n1() && v1.n2() == v2.n2();
	}
	else if constexpr (ndims<T1>() == 3 && ndims<T2>() == 3) {
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
	for (Long i = 0; i < v1.n1(); ++i)
		for (Long j = 0; j < v1.n2(); ++j)
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

// copy one row of dense matrix to a dense vector
template <class Tvec, class Tmat,
	SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_row(Tvec &v, const Tmat &a, Long_I row)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_SHAPE
	if (v.size() != Nc) {
		SLS_ERR("wrong shape!");
	}
#endif
	if constexpr (is_rmajor<Tmat>()) { // row major
		veccpy(v.ptr(), a.ptr() + Nc*row, Nc);
	}
	else if constexpr (is_cmajor<Tmat>()) { // column major
		auto p = a.ptr() + row;
		for (Long i = 0; i < Nc; ++i) {
			v[i] = *p;
			p += Nr;
		}
	}
	else
		SLS_ERR("unknown!");
}

// copy a dense vector to one row of dense matrix
template <class Tvec, class Tmat,
	SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
	inline void copy_row(Tmat &a, const Tvec &v, Long_I row)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_SHAPE
	if (v.size() != Nc) {
		SLS_ERR("wrong shape!");
	}
#endif
	if constexpr (is_rmajor<Tmat>()) { // row major
		veccpy(a.ptr() + Nc * row, v.ptr(), Nc);
	}
	else if constexpr (is_cmajor<Tmat>()) { // column major
		auto p = a.ptr() + row;
		for (Long i = 0; i < Nc; ++i) {
			*p = v[i];
			p += Nr;
		}
	}
	else
		SLS_ERR("unknown!");
}

// copy one column of dense matrix to a dense vector
template <class Tvec, class Tmat,
	SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_col(Tvec &v, const Tmat &a, Long_I col)
{
#ifdef SLS_CHECK_SHAPE
	if (v.size() != a.n1()) {
		SLS_ERR("wrong shape!");
	}
#endif
	Long Nr = a.n1(), Nc = a.n2();
	if constexpr (is_cmajor<Tmat>()) { // column major
		veccpy(v.ptr(), a.ptr() + Nr*col, Nr);
	}
	else if constexpr (is_rmajor<Tmat>()) { // row major
		auto p = a.ptr() + col;
		for (Long i = 0; i < Nr; ++i) {
			v[i] = *p;
			p += Nc;
		}
	}
	else
		SLS_ERR("unknown!");
}

// copy a dense vector to one column of dense matrix
template <class Tvec, class Tmat,
	SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_col(Tmat &a, const Tvec &v, Long_I col)
{
#ifdef SLS_CHECK_SHAPE
	if (v.size() != a.n1()) {
		SLS_ERR("wrong shape!");
	}
#endif
	Long Nr = a.n1(), Nc = a.n2();
	if constexpr (is_cmajor<Tmat>()) { // column major
		veccpy(a.ptr() + Nr * col, v.ptr(), Nr);
	}
	else if constexpr (is_rmajor<Tmat>()) { // row major
		auto p = a.ptr() + col;
		for (Long i = 0; i < Nr; ++i) {
			*p = v[i];
			p += Nc;
		}
	}
	else
		SLS_ERR("unknown!");
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

template <class Tv, class T1, class T2, SLS_IF(is_dense<Tv>())>
inline void linspace(Tv &v, const T1 &first, const T2 &last)
{
	linspace_vss(v.ptr(), (contain_type<Tv>)first, (contain_type<Tv>)last, v.size());
}

// === vectorized math functions ===

template <class T, SLS_IF(is_dense<T>())>
void sqrt(T &v)
{
	sqrt_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T,T1>())>
void sqrt(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	sqrt_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class Ts, SLS_IF(
	is_dense<T>() && is_scalar<Ts>())>
void pow(T &v, const Ts &s)
{
	pow_vs(v.ptr(), s, v.size());
}

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T, T1>() && is_scalar<T2>())>
void pow(T &v, const T1 &v1, const T2 &s)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	pow_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void sin(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	sin_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void cos(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	cos_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, SLS_IF(is_fpt_dense<T>())>
void exp(T &v)
{
	exp_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void exp(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	exp_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
void tan(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	tan_vv(v.ptr(), v1.ptr(), v1.size());
}

// === vector/matrix arithmetics ===

// flip(v)
template <class T, SLS_IF(is_dense_vec<T>())>
inline void flip(T &v)
{
	flip_v(v.ptr(), v.size());
}

// v = flip(v)
template <class T, class T1, SLS_IF(
	is_dense_vec<T>() && is_dense_vec<T1>()
)>
inline void flip(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	flip(v.ptr(), v1.ptr(), v1.size());
}

// matrix transpose

// trans(v)
template <class T, SLS_IF(is_dense_mat<T>())>
inline void trans(T &v)
{
#ifdef SLS_CHECK_SHAPE
	if (v.n1() != v.n2())
		SLS_ERR("illegal shape!");
#endif
	for (Long i = 0; i < v.n1(); ++i)
		for (Long j = 0; j < i; ++j)
			swap(v(i, j), v(j, i));
}

// v = trans(v)
// TODO: this is inefficient
template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>()
)>
inline void trans(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (v.n1() != v1.n2() || v.n2() != v1.n1())
		SLS_ERR("wrong size!");
#endif
	for (Long i = 0; i < v.n1(); ++i)
		for (Long j = 0; j < v.n2(); ++j)
			v(i, j) = v1(j, i);
}

// hermitian conjugate

// her(v)
template <class T, SLS_IF(is_comp_dense<T>() && is_dense_mat<T>())>
inline void her(T &v)
{
#ifdef SLS_CHECK_SHAPE
	if (v.n1() != v.n2()) SLS_ERR("illegal shape!");
#endif
	// diagonal
	for (Long i = 0; i < v.n1(); ++i)
		v(i, i) = conj(v(i, i));
	// off-diagonal
	for (Long i = 0; i < v.n1(); ++i) {
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
	is_comp_dense<T>() && is_comp_dense<T1>() &&
	type_num<contain_type<T>>() >= type_num<contain_type<T1>>()
)>
inline void her(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (v.n1() != v1.n2() || v.n2() != v1.n1())
		SLS_ERR("wrong shape!");
#endif
	for (Long i = 0; i < v.n1(); ++i) {
		for (Long j = 0; j < v.n2(); ++j)
			v(i, j) = conj(v1(j, i));
	}
}

// v += v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator+=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	plus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v -= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator-=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	minus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v *= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator*=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	times_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v /= v

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void operator/=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
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

template <class T, class Ts, SLS_IF(is_Dvector<T>() && is_scalar<Ts>())>
inline void operator/=(T &v, const Ts &s)
{
	divide_equals_vs(v.ptr(), s, v.size(), v.step());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Plus(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v + v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Plus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2) || !shape_cmp(v, v1)) SLS_ERR("wrong shape!");
#endif
	plus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	minus_vv(v.ptr(), v1.ptr(), v1.size());
}

// v = s - v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Minus(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	minus_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v - s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Minus(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	minus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v - v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Minus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2) || !shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	minus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

// v = v * s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Times(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s * v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Times(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v * v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Times(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2) || !shape_cmp(v, v1)) SLS_ERR("wrong shape!");
#endif
	times_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// v = v / s

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Divide(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	divide_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s / v

template <class T, class T1, class Ts, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_scalar<Ts>()
)>
inline void Divide(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	divide_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v / v

template <class T, class T1, class T2, SLS_IF(
	is_dense<T>() && is_same_contain<T,T1>() && is_same_contain<T,T2>()
)>
inline void Divide(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2) || !shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	divide_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	real_vv(v.ptr(), v1.ptr(), v1.size());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	imag_vv(v.ptr(), v1.ptr(), v1.size());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	abs_vv(v.ptr(), v1.ptr(), v1.size());
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
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong size!");
#endif
	conj_vv(v.ptr(), v1.ptr(), v1.size());
}

// dot products ( sum conj(v1[i])*v2[i] )

// s = dot(v, v)
template <class T1, class T2, SLS_IF(is_dense_vec<T1>() && is_dense_vec<T2>())>
inline auto dot(const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v1, v2))
		SLS_ERR("wrong shape!");
#endif
	return dot_vv(v1.ptr(), v2.ptr(), v2.size());
}

// outer product ( v(i,j) = conj(v1[i])*v2[j] )
// outprod(v, v, v)
template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_dense_vec<T1>() && is_dense_vec<T2>()
)>
inline void outprod(T &v, const T1 &v1, const T2 &v2)
{
	SLS_ERR("TODO");
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

template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_dense_vec<T1>() && is_dense_vec<T2>()
)>
inline void outprod_par(T &v, const T1 &v1, const T2 &v2)
{
	SLS_ERR("TODO");
	/*Long N1{ v1.size() }, N2{ v2.size() };
	v.resize(N1, N2);
#pragma omp parallel for
	for (Long i = 0; i < N1; ++i) {
		Comp *pc, v1_i;
		pc = v.ptr(i);
		v1_i = v1[i];
		for (Long j = 0; j < N2; ++j)
			pc[j] = v1_i*v2[j];
	}*/
}

// matrix-vector multiplication
template <class T, class T1, class T2, SLS_IF(
	is_dense_vec<T>() && is_dense_mat<T1>() && is_dense_vec<T2>()
)>
inline void mul(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != x.size() || y.size() != a.n1())
		SLS_ERR("illegal shape!");
#endif
	Long i, j, Nr_a = a.n1(), Nc_a = a.n2();
	vecset(y.ptr(), contain_type<T>(), Nr_a);
	for (i = 0; i < Nr_a; ++i) {
		for (j = 0; j < Nc_a; ++j)
			y[i] += a(i, j) * x[j];
	}
}

// parallel version
template <class T, class T1, class T2, SLS_IF(
	is_dense_vec<T>() && is_dense_mat<T1>() && is_dense_vec<T2>()
)>
inline void mul_par(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != x.size() || y.size() != a.n1())
		SLS_ERR("illegal shape!");
#endif
	Long Nr_a = a.n1(), Nc_a = a.n2();
	vecset(y.ptr(), contain_type<T>(), Nr_a);
#pragma omp parallel for
	for (Long i = 0; i < Nr_a; ++i) {
		for (Long j = 0; j < Nc_a; ++j)
			y[i] += a(i, j) * x[j];
	}
}

// vector-matrix multiplication (row vector assumed)
template <class T, class T1, class T2, SLS_IF(
	is_dense_vec<T>() && is_dense_vec<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &x, const T2 &a)
{
	Long Nr_a = a.n1(), Nc_a = a.n2();
#ifdef SLS_CHECK_SHAPE
	if (x.size() != a.n1() || y.size() != Nc_a)
		SLS_ERR("illegal shape!");
#endif
	vecset(y.ptr(), contain_type<T>(), Nc_a);
	for (Long j = 0; j < Nc_a; ++j) {
		for (Long i = 0; i < Nr_a; ++i)
			y[j] += x[i] * a(i, j);
	}
}

// parallel version
template <class T, class T1, class T2, SLS_IF(
	is_dense_vec<T>() && is_dense_vec<T1>() && is_dense_mat<T2>()
)>
inline void mul_par(T &y, const T1 &x, const T2 &a)
{
	Long Nr_a = a.n1(), Nc_a = a.n2();
#ifdef SLS_CHECK_SHAPE
	if (x.size() != a.n1() || y.size() != Nc_a)
		SLS_ERR("illegal shape!");
#endif
	vecset(y.ptr(), contain_type<T>(), Nc_a);
#pragma omp parallel for
	for (Long j = 0; j < Nc_a; ++j) {
		for (Long i = 0; i < Nr_a; ++i)
			y[j] += x[i] * a(i, j);
	}
}

// matrix-matrix multiplication

template <class T, class T1, class T2, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &a, const T2 &x)
{
	Long Nr_a = a.n1(), Nc_a = a.n2(), Nc_x = x.n2();
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != x.n1() || y.n1() != Nr_a || y.n2() != Nc_x)
		SLS_ERR("illegal shape!");
#endif
	vecset(y.ptr(), contain_type<T>(), Nr_a*Nc_x);
	for (Long i = 0; i < Nr_a; ++i) {
		for (Long j = 0; j < Nc_x; ++j) {
			for (Long k = 0; k < Nc_a; ++k)
				y(i, j) += a(i, k) * x(k, j);
		}
	}
}

// === numerical integration ===

// indefinite integral;
// use cumsum(y)*dx instead
template <class T, class T1, SLS_IF(is_dense_vec<T>() && is_dense_vec<T1>())>
void cumsum(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("illegal shape!");
#endif
		cumsum_vv(v.ptr(), v1.ptr(), v1.size());
}

} // namespace slisc
