// low-level arithmetic
// use pointers for array input/output

#include "scalar_arith.h"

namespace slisc {

// array comparison
template <class T1, class T2>
Bool equals_to_vv(const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		if (v1[i] != v2[i])
			return false;
	return true;
}

template <class T1, class T2>
Bool equals_to_vs(const T1 *v, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		if (v[i] != s)
			return false;
	return true;
}

// +=
template <class T, class T1>
inline void plus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

template <class T, class T1>
inline void plus_equals_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += s;
}

// -=
template <class T, class T1>
inline void minus_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

template <class T, class T1>
inline void minus_equals_vs(T *v, const T1 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= s;
}

// *=
template <class T, class T1>
inline void times_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

template <class T, class T1>
inline void times_equals_vs(T *v, const T1 &s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

// v /= v
template <class T, class T1>
inline void divide_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= v1[i];
}

template <class T, class T1>
inline void divide_equals_vs(T *v, const T1 &s, Long_I N)
{
	if constexpr (std::is_floating_point<T>::value) {
		times_equals_vs(v, 1./s, N);
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] /= s;
	}
}

// operator%
template <class T, class T1, class T2>
inline void rem_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] % s;
}

template <class T, class T1, class T2>
inline void mod_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = mod(v1[i], s);
}

template <class T, class T1, class T2>
inline void plus_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + s;
}

template <class T, class T1, class T2>
inline void plus_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

// minus

template <class T, class T1>
inline void minus_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v1[i];
}

template <class T, class T1, class T2>
inline void minus_vsv(T *v, const T1 &s, const T2 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

template <class T, class T1, class T2>
inline void minus_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

template <class T, class T1, class T2>
inline void minus_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - v2[i];
}

template <class T, class T1, class T2>
void times_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i) {
		v[i] = v1[i] * s;
	}
}

template <class T, class T1, class T2>
inline void times_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] * v2[i];
}

template <class T, class T1, class T2>
inline void divide_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	T2 sInv{ 1 / s };
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] * sInv;
}

template <class T, class T1, class T2>
inline void divide_vsv(T *v, const T1 &s, const T2 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s / v1[i];
}

template <class T, class T1, class T2>
inline void divide_vvv(T *v, const T1 *v1, const T2 *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] / v2[i];
}

template <class T, class Tc>
inline void real_vv(T *v, const std::complex<Tc> *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = real(v1[i]); 
}

template <class T, class Tc>
inline void imag_vv(T *v, const std::complex<Tc> *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = imag(v1[i]);
}

template <class T, class T1>
inline void abs_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = abs(v1[i]);
}

template <class Tc, class T>
inline void to_comp_vv(std::complex<Tc> *v, const T *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i];
}

template <class T1, class T2, class T = typename promo_type<T1, T2>::type>
inline void dot_svv(T &s, const T1 *v1, const T2 *v2, Long_I N)
{
	s = T();
	for (Long i = 0; i < N; ++i) {
		if constexpr (is_complex<T1>())
			s += std::conj(v1[i]) * v2[i];
		else
			s += v1[i] * v2[i];
	}
}

// sqrt(v)
template <class T, class T1>
inline void sqrt_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		for (Long i = 0; i < N; ++i)
			v[i] = sqrt(T(v1[i]));
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] = sqrt(v1[i]);
	}
}

template <class T, class T1>
inline void invSqrt_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		// if T is the smaller type
		// convert v1[i] to T before sqrt()
		typedef typename rm_complex<T>::type Tr;
		constexpr Tr one = Tr(1);
		for (Long i = 0; i < N; ++i)
			v[i] = one / sqrt(T(v1[i]));
	}
	else {
		// if T is the larger type
		typedef typename rm_complex<T1>::type Tr;
		constexpr Tr one = Tr(1);
		for (Long i = 0; i < N; ++i)
			v[i] = one / sqrt(v1[i]);
	}
}

template <class T, class T1>
inline void sin_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		for (Long i = 0; i < N; ++i)
			v[i] = sin(T(v1[i]));
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] = sin(v1[i]);
	}
}

template <class T, class T1>
inline void cos_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		for (Long i = 0; i < N; ++i)
			v[i] = cos(T(v1[i]));
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] = cos(v1[i]);
	}
}

template <class T, class T1>
inline void exp_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		for (Long i = 0; i < N; ++i)
			v[i] = exp(T(v1[i]));
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] = exp(v1[i]);
	}
}

template <class T, class T1>
inline void tan_vv(T *v, const T1 *v1, Long_I N)
{
	if constexpr (type_num<T>() < type_num<T1>()) {
		for (Long i = 0; i < N; ++i)
			v[i] = tan(T(v1[i]));
	}
	else {
		for (Long i = 0; i < N; ++i)
			v[i] = tan(v1[i]);
	}
}

} // nemaspace slisc
