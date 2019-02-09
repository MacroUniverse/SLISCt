// low-level arithmetic
// use pointers for array input/output

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

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>()
)>
inline void plus_equals_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += s;
}

// v += v

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>() ||
	is_Comp<T>() && is_Doub<T1>()
)>
inline void plus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

// v -= s

template <class T, class T1, SLS_IF(
	is_Char<T>() && is_Char<T1>() ||
	is_Int<T>() && is_Int<T1>() ||
	is_Llong<T>() && is_Llong<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>() ||
	is_Comp<T>() && is_Doub<T1>()
)>
inline void minus_equals_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= s;
}

// v -= v

template <class T, class T1, SLS_IF(
	is_Char<T>() && is_Char<T1>() ||
	is_Int<T>() && is_Int<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>() ||
	is_Comp<T>() && is_Doub<T1>()
)>
inline void minus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

// v *= s

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>()
)>
inline void times_equals_vs(T *v, const T1 &s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

// v *= v

template <class T, class T1, SLS_IF(
	is_Char<T>() && is_Char<T1>() ||
	is_Int<T>() && is_Int<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>() ||
	is_Comp<T>() && is_Doub<T1>()
)>
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
	is_Comp<T>() && is_Doub<T1>()
)>
inline void divide_equals_vs(T *v, const T1 &s, Long_I N)
{ times_equals_vs(v, INV(s), N); }

// v /= v

template <class T, class T1, SLS_IF(
	is_Int<T>() && is_Int<T1>() ||
	is_Float<T>() && is_Float<T1>() ||
	is_Doub<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Doub<T1>() ||
	is_Comp<T>() && is_Comp<T1>()
)>
inline void divide_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= v1[i];
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
inline void rem_vs(T *v, const T &s, Long_I N)
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
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void plus_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

// -v
template <class T, SLS_IF(
	is_Int<T>() || is_Float<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void minus_v(T *v, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= -1;
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
	is_Float<T>() && is_Float<T1>() && is_Float<T2>() ||
	is_Doub<T>() && is_Doub<T1>() && is_Doub<T2>() ||
	is_Comp<T>() && is_Doub<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>() ||
	is_Comp<T>() && is_Comp<T1>() && is_Comp<T2>()
)>
inline void divide_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	times_vvs(v, v1, INV(s), N);
}

// v = s / v
template <class T, class T1, class T2, SLS_IF(
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

inline void real_v(Comp *v, Long_I N)
{
	Doub *pd = (Doub *)v;
	for (Long i = 1; i < N; i += 2)
		pd[i] = 0.;
}

// v = real(v)

inline void real_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = real(v1[i]); 
}

// imag(v)

inline void imag_v(Comp *v, Long_I N)
{
	Doub *pd = (Doub *)v;
	for (Long i = 0; i < N; i += 2)
		pd[i] = 0.;
}

// v = imag(v)

inline void imag_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = imag(v1[i]);
}

// abs(v)

template <class T, SLS_IF(
	is_Int<T>() && is_Doub<T>() && is_Comp<T>()
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

// v = comp(v)

template <class T, class T1, SLS_IF(
	is_Comp<T>() && is_Float<T1>() ||
	is_Comp<T>() && is_Doub<T1>()
)>
inline void to_comp_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i];
}

// s = sum(v)

template <class T, class T1, SLS_IF(
	is_Llong<Long>() &&
	(is_Bool<T>() || is_Int<T>() || is_Llong<T>())
)>
inline Long sum_v(const T *v, Long_I N)
{
	Long s = v[0];
	for (Long i = 1; i < N; ++i)
		s += v[i];
	return s;
}

template <class T, SLS_IF(
	is_Float<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline T sum_v(const T *v, Long_I N)
{
	T s = v[0];
	for (Long i = 1; i < N; ++i)
		s += v[i];
	return s;
}

// s = max(v)

template <class T, SLS_IF(
	is_Int<T>() || is_Float<T>() ||is_Doub<T>()
)>
inline T max_v(const T *v, Long_I N)
{
	T s = v[0], val;
	for (Long i = 1; i < N; ++i) {
		if (s < v[i])
			s = v[i];
	}
	return s;
}

// s = max_abs(v)

template <class T, SLS_IF( is_real<T>() )>
inline T max_abs(const T *v, Long_I N)
{
	T s = abs(v[0]), val;
	for (Long i = 1; i < N; ++i) {
		val = abs(v[i]);
		if (s < val)
			s = val;
	}
	return s;
}

template <class Tc, SLS_IF(is_comp<Tc>())>
inline rm_comp<Tc> max_abs_v(const Tc *v, Long_I N)
{
	rm_comp<Tc> s = abs(v[0]), val;
	for (Long i = 1; i < N; ++i) {
		val = abs(v[i]);
		if (s < val)
			s = val;
	}
	return s;
}

// conj(v)

inline void conj(Comp *v, Long_I N)
{
	Doub *p = (Doub *)v;
	for (Long i = 1; i < N; i += 2)
		p[i] = -p[i];
}

// s = dot(v, v)

inline Doub dot_vv(const Doub *v1, const Char *v2, Long_I N)
{
	Doub s = v1[0] * v2[0];
	for (Long i = 1; i < N; ++i) {
		s += v1[i] * v2[i];
	}
	return s;
}

inline Doub dot_vv(const Doub *v1, const Doub *v2, Long_I N)
{
	Doub s = v1[0] * v2[0];
	for (Long i = 1; i < N; ++i) {
		s += v1[i] * v2[i];
	}
	return s;
}

inline Comp dot_vv(const Doub *v1, const Comp *v2, Long_I N)
{
	Comp s = v1[0] * v2[0];
	for (Long i = 1; i < N; ++i) {
		s += v1[i] * v2[i];
	}
	return s;
}

inline Comp dot_vv(const Comp *v1, const Doub *v2, Long_I N)
{
	Comp s = conj(v1[0]) * v2[0];
	for (Long i = 1; i < N; ++i) {
		s += conj(v1[i]) * v2[i];
	}
	return s;
}

inline Comp dot_vv(const Comp *v1, const Comp *v2, Long_I N)
{
	Comp s = conj(v1[0]) * v2[0];
	for (Long i = 1; i < N; ++i) {
		s += conj(v1[i]) * v2[i];
	}
	return s;
}

template <class T>
inline void flip(T *v, Long_I N)
{
	for (Long i = 0; i < N / 2; ++i)
		swap(v[i], v[N - i - 1]);
}

template <class T, class T1, SLS_IF(is_scalar<T>() && is_scalar<T1>())>
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

inline void sqrt_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sqrt(v1[i]);
}

inline void sqrt_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sqrt(v1[i]);
}

// v = 1/sqrt(v)

inline void invSqrt_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = pow(v1[i], -0.5f);
}

inline void invSqrt_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = pow(v1[i], -0.5f);
}

// v = sin(v)

inline void sin_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sin(v1[i]);
}

inline void sin_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = sin(v1[i]);
}

// v = cos(v)

inline void cos_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = cos(v1[i]);
}

inline void cos_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = cos(v1[i]);
}

// v = exp(v)

inline void exp_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = exp(v1[i]);
}

inline void exp_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = exp(v1[i]);
}

// v = tan(v)

inline void tan_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = tan(v1[i]);
}

inline void tan_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = tan(v1[i]);
}

// v = cumsum(v)
inline void cumsum_vv(Doub *v, const Doub *v1, Long_I N)
{
	v[0] = v1[0];
	for (Long i = 1; i < N; ++i)
		v[i] = v[i - 1] + v1[i];
}

} // nemaspace slisc
