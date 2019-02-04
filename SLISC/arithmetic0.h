// low-level arithmetic
// use pointers for array input/output

#include "slisc.h"

namespace slisc {

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
inline void divide_equals_vv(T &v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= v1[i];
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
inline void times_vvs(T *v, const T1 *v1, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v(i) = v1(i) * s;
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
inline void divide_vvv(T *v, const T1 *v1, const T2 *v2)
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

inline void doub2comp0_vv(Comp *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i];
}

template <class T, class T1, class T2>
inline const T dot_vv(const T1 *v1, const T2 *v2, Long_I N)
{
	T s();
	for (Long i = 0; i < N; ++i)
		// std::conj() does not have overhead for non-complex types
		s += std::conj(v1[i]) * v2[i];
	return s;
}

// vector-matrix multiplication
template <class T, class T1, class T2>
inline void mul_vvm(T *y, const T1 *x, const T2 *a, Long_I a_Nr, Long_I a_Nc)
{
	Long i, j;
	vecset(y, T(), a_Nr);
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i)
			y[j] += x[i] * a[rsub2ind(a_Nc, i, j)];
	}
}

// Cmat-vector multiplications (slow)
template <class T, class T1, class T2>
inline void cmul_vmv(T *y, const T1 *a, Long_I a_Nr, Long_I a_Nc, const T2 *x)
{
	Long i, j;
	vecset(y, T(), a_Nr);
	for (i = 0; i < a_Nr; ++i) {
		for (j = 0; j < a_Nr; ++j)
			y[i] += a[csub2ind(a_Nr, i, j)] * x[j];
	}
}

// Cmat-Cmat multiplication (slow)
template <class T, class T1, class T2> // column major
inline void cmul_mmm(T *c, const T1 *a, Long_I Nr_a, Long_I Nc_a, const T2 *b, Long_I Nc_b)
{
	Long i, j, k;
	vecset(c, T(), Nr_a);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			for (k = 0; k < Nk; ++k)
				c[csub2ind(Nr_a, i, j)] += a[csub2ind(Nr_a, i, j)] * b[csub2ind(Nc_a, k, j)];
		}
	}
}

// sqrt(v)
template <class T, class T1>
inline void sqrt_vv(Vbase<T> &v, const Vbase<T1> &v1, Long_I N)
{
	Long i, N{ v1.size() };
	for (i = 0; i < N; ++i)
		v[i] = std::sqrt(v1[i]);
}

template <class T, class T1>
inline void invSqrt_vv(T *v, const T1 *v1, Long_I N)
{
	for (i = 0; i < N; ++i)
		v[i] = 1./std::sqrt(v1[i]);
}

template <class T, class T1>
inline void sin_vv(Vbase<T> &v, const Vbase<T1> &v1, Long_I N)
{
	for (i = 0; i < N; ++i)
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
