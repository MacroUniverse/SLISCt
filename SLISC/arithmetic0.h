// low-level arithmetic
// use pointers for array input/output

#include "slisc.h"
#include "scalar_arith.h"

namespace slisc {

// array copying
template<class T>
inline void vecset(T *dest, const T &val, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = val;
}

template<class T>
inline void veccpy(T *dest, const T *src, Long_I n)
{
	memcpy(dest, src, n * sizeof(T));
}

template<class T, class T1>
inline void veccpy(T *dest, const T1 *src, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = src[i];
}

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
template <class T>
inline void rem_vvs(T *v, const T *v1, const T &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] % s;
}

// mod(v, v, s)
Int mod(Int_I i, Int_I n);
Long mod(Long_I i, Long_I n);

template <class T>
inline void mod_vvs(T *v, const T *v1, const T &s, Long_I N)
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

inline void real_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::real(v1[i]);
}

inline void imag_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::imag(v1[i]);
}

template <class T, class T1>
inline void abs_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = abs(v1[i]);
}

inline void doub2comp_vv(Comp *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i];
}

template <class T, class T1, class T2>
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
inline void sqrt_vv(Vbase<T> &v, const Vbase<T1> &v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::sqrt(v1[i]);
}

template <class T, class T1>
inline void invSqrt_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = 1./std::sqrt(v1[i]);
}

template <class T, class T1>
inline void sin_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::sin(v1[i]);
}

template <class T, class T1>
inline void cos_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::cos(v1[i]);
}

template <class T, class T1>
inline void exp_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::exp(v1[i]);
}

template <class T, class T1>
inline void tan_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = std::tan(v1[i]);
}

} // nemaspace slisc
