// low-level arithmetic
// use pointers for array input/output
#pragma once
#include "scalar_arith.h"

namespace slisc {

// v == s

template <class T1, class T2>
Bool equals_to_vs(const T1 *v, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		if (v[i] != s)
			return false;
	return true;
}

// v == v

template <class T1, class T2>
Bool equals_to_vv(const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		if (v1[i] != v2[i])
			return false;
	return true;
}

// v += s

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void plus_equals_vs(T *v, const T1 &s, Long_I N)
{
	T s1 = T(s);
	for (Long i = 0; i < N; ++i)
		v[i] += s1;
}

// v += v

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void plus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

// v -= s

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void minus_equals_vs(T *v, const T1 &s, Long_I N)
{
	T s1 = (T)s;
	for (Long i = 0; i < N; ++i)
		v[i] -= s1;
}

// v -= v

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void minus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

// v *= s

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void times_equals_vs(T *v, const T1 &s, Long N)
{
	T s1 = T(s);
	for (Long i = 0; i < N; ++i)
		v[i] *= s1;
}

// v *= v

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void times_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

// v /= s

template <class T, class T1, SLS_IF(
	is_Char<T>() && is_Char<T1>() ||
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>()
)>
inline void divide_equals_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= s;
}

template <class T, class T1, SLS_IF(
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>()
)>
inline void divide_equals_vs(T *v, const T1 &s, Long_I N)
{
	times_equals_vs(v, INV(s), N);
}

// v /= v

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void divide_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= v1[i];
}

// mod(v, s)
template <class T, class T1, SLS_IF(
	is_Char<T>() && is_Char<T1>() ||
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>()
)>
inline void mod_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = mod(v[i], s);
}

// v = mod(v, s)
template <class T, class T1, class T2, SLS_IF(
	is_Char<T>() && is_Char<T1>() && is_Char<T2>() ||
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Llong<T>() && is_Llong<T1>() && is_Llong<T2>()
)>
inline void mod_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = mod(v1[i], s);
}

// v %= s

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>()
)>
inline void rem_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] %= s;
}

// v = v % s

template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Llong<T>() && is_Llong<T1>() && is_Llong<T2>()
)>
inline void rem_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] % s;
}

// v = v + s

template <class T, class T1, class T2, SLS_IF(
	is_Char<T>() && is_Char<T1>() && is_Char<T2>() ||
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void plus_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + s;
}

// v = v + v

template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Llong<T>() && is_Llong<T1>() && is_Llong<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>() ||
	is_Lcomp<T>() && is_Lcomp<T1>() && is_Lcomp<T2>()
)>
inline void plus_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

// -v
template <class T, SLS_IF(
	is_Int<T>() || is_Llong<T>() || is_Float<T>() ||
	is_Doub<T>() || is_Comp<T>()
)>
inline void minus_v(T *v, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v[i];
}

// v = -v

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>()
)>
inline void minus_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v1[i];
}

// v = s - v

template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void minus_vsv(T *v, const T1 &s, const T2 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

// v = v - s
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void minus_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

// v = v - v
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void minus_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - v2[i];
}

// v = v * s
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
void times_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i) {
		v[i] = v1[i] * s;
	}
}

// v = v * v
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void times_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] * v2[i];
}

// v = v / s

template <class T, class T1, class T2, SLS_IF(
	is_Char<T>() && is_Char<T1>() && is_Char<T2>() ||
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Llong<T>() && is_Llong<T1>() && is_Llong<T2>()
)>
inline void divide_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] / s;
}

template <class T, class T1, class T2, SLS_IF(
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void divide_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	times_vvs(v, v1, INV(s), N);
}

// v = s / v
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Llong<T>() && is_Llong<T1>() && is_Llong<T2>() ||
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()

)>
inline void divide_vsv(T *v, const T1 &s, const T2 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s / v1[i];
}

// v = v / v
template <class T, class T1, class T2, SLS_IF(
	is_Int<T>() && is_Int<T1>() && is_Int<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void divide_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] / v2[i];
}

// real(v)

template <class T, SLS_IF(is_comp<T>())>
inline void real_v(T *v, Long_I N)
{
	rm_comp<T> *pr = (rm_comp<T> *)v;
	for (Long i = 1; i < 2*N; i += 2)
		pr[i] = 0;
}

// v = real(v)

template <class Tr, class Tc, SLS_IF(is_comp<Tc>() && is_same<complex<Tr>, Tc>())>
inline void real_vv(Tr *v, const Tc *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = real(v1[i]); 
}

// imag(v)

template <class T, SLS_IF(is_comp<T>())>
inline void imag_v(T *v, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = imag(v[i]);
}

// v = imag(v)

template <class Tr, class Tc, SLS_IF(is_comp<Tc>() && is_same<complex<Tr>, Tc>())>
inline void imag_vv(Tr *v, const Tc *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = imag(v1[i]);
}

// abs(v)

template <class T, SLS_IF(
	is_Int<T>() || is_Float<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void abs_v(T *v, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = abs(v[i]);
}

// v = abs(v)

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Doub<T>() && is_Comp<T1>()
)>
inline void abs_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = abs(v1[i]);
}

// s = sum(v)

template <class T, SLS_IF(
	is_Llong<Long>() && is_integral<T>()
)>
inline Long sum_v(const T *v, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	Long s = v[0];
	for (Long i = 1; i < N; ++i)
		s += v[i];
	return s;
}

template <class T, SLS_IF(
	is_floating_point<T>() || is_comp<T>()
)>
inline T sum_v(const T *v, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	T s = v[0];
	for (Long i = 1; i < N; ++i)
		s += v[i];
	return s;
}

// s = max(v)

template <class T, SLS_IF(is_real<T>() && !is_Bool<T>())>
inline T max_v(const T *v, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	T s = v[0], val;
	for (Long i = 1; i < N; ++i) {
		if (s < v[i])
			s = v[i];
	}
	return s;
}

// s = max_abs(v)

template <class T, SLS_IF(is_scalar<T>() && !is_Bool<T>())>
inline rm_comp<T> max_abs_v(const T *v, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	rm_comp<T> s = abs(v[0]), val;
	for (Long i = 1; i < N; ++i) {
		val = abs(v[i]);
		if (s < val)
			s = val;
	}
	return s;
}

// conj(v)

template <class T, SLS_IF(is_comp<T>())>
inline void conj_v(T *v, Long_I N)
{
	rm_comp<T> *p = (rm_comp<T> *)v;
	for (Long i = 1; i < 2*N; i += 2)
		p[i] = -p[i];
}

// v = conj(v)

template <class T, class T1, SLS_IF(
	is_comp<T1>() && is_promo<T, T1>()
)>
inline void conj_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = conj(v1[i]);
}

// s = dot(v, v)
template <class T1, class T2, SLS_IF(
	is_scalar<T1>() && is_scalar<T2>()
)>
inline auto dot_vv(const T1 *v1, const T2 *v2, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	auto s = conj(v1[0]) * v2[0];
	for (Long i = 1; i < N; ++i) {
			s += conj(v1[i]) * v2[i];
	}
	return s;
}

template <class T>
inline void flip_v(T *v, Long_I N)
{
	for (Long i = 0; i < N / 2; ++i)
		swap(v[i], v[N - i - 1]);
}

template <class T, class T1, SLS_IF(
	is_same<T, T1>() ||
	!is_same<T, T1>() && is_scalar<T>() && is_scalar<T1>())>
inline void flip(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[N - i - 1];
}

template <class T, SLS_IF(is_scalar<T>())>
inline void linspace_vss(T *v, const T &first, const T &last, Long N)
{
	T delta = (last - first) / (N - 1);
	for (Long i = 0; i < N; ++i)
		v[i] = first + delta * i;
}

// v = sqrt(v)
template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void sqrt_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sqrt(v1[i]);
}

// v = 1/sqrt(v)

template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void invSqrt_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = pow(v1[i], rm_comp<T>(-0.5));
}

// v = sin(v)

template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void sin_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sin(v1[i]);
}

// v = cos(v)

template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void cos_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = cos(v1[i]);
}

// v = exp(v)

template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void exp_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = exp(v1[i]);
}

// v = tan(v)

template <class T, SLS_IF(is_floating_point<T>() || is_comp<T>())>
inline void tan_vv(T *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = tan(v1[i]);
}

// v = cumsum(v)

template <class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void cumsum_vv(T *v, const T1 *v1, Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
	if (N <= 0) error("illegal length!");
#endif
	v[0] = v1[0];
	for (Long i = 1; i < N; ++i)
		v[i] = v[i - 1] + v1[i];
}

} // nemaspace slisc
