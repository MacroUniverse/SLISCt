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

// v += v
inline void plus_equals_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

inline void plus_equals_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

inline void plus_equals_vv(Comp *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

inline void plus_equals_vv(Comp *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += v1[i];
}

// v += s
inline void plus_equals_vs(Float *v, Float_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += s;
}

inline void plus_equals_vs(Doub *v, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += s;
}

inline void plus_equals_vs(Comp *v, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] += s;
}

// v -= v

inline void minus_equals_vv(Float *v, Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

inline void minus_equals_vv(Doub *v, Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

inline void minus_equals_vv(Comp *v, Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

inline void minus_equals_vv(Comp *v, Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= v1[i];
}

// v -= s

inline void minus_equals_vs(Doub *v, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= s;
}

inline void minus_equals_vs(Comp *v, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= s;
}

inline void minus_equals_vs(Comp *v, Comp_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] -= s;
}

// v *= v
inline void times_equals_vv(Float *v, Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

inline void times_equals_vv(Doub *v, Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

inline void times_equals_vv(Comp *v, Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

inline void times_equals_vv(Comp *v, Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= v1[i];
}

// v *= s
inline void times_equals_vs(Float *v, Float_I s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

inline void times_equals_vs(Doub *v, Doub_I s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

inline void times_equals_vs(Comp *v, Doub_I s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

inline void times_equals_vs(Comp *v, Comp_I s, Long N)
{
	for (Long i = 0; i < N; ++i)
		v[i] *= s;
}

// v /= v
inline void divide_equals_vv(T *v, const T1 *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= v1[i];
}

// v /= s
inline void divide_equals_vs(Float *v, const Float &s, Long_I N)
{ times_equals_vs(v, 1.f/s, N); }

inline void divide_equals_vs(Doub *v, const Doub &s, Long_I N)
{ times_equals_vs(v, 1./s, N); }

inline void divide_equals_vs(Comp *v, const Doub &s, Long_I N)
{ times_equals_vs(v, 1./s, N); }

inline void divide_equals_vs(Char *v, const Char &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= s;
}

inline void divide_equals_vs(Int *v, const Int &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= s;
}

inline void divide_equals_vs(Llong *v, const Llong &s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] /= s;
}

// operator%

inline void mod_vvs(Char *v, const Char *v1, Char_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = mod(v1[i], s);
}

inline void rem_vvs(Int *v, const Int *v1, Int_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] % s;
}

inline void rem_vvs(Long *v, const Long *v1, Long_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] % s;
}

inline void plus_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + s;
}

inline void plus_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] + v2[i];
}

// minus

inline void minus_vv(Int *v, const Int *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v1[i];
}

inline void minus_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v1[i];
}

inline void minus_vv(Comp *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = -v1[i];
}

inline void minus_vsv(Int *v, Int_I s, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

inline void minus_vsv(Doub *v, Doub_I s, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

inline void minus_vsv(Comp *v, Comp_I s, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s - v1[i];
}

// v = v - s

inline void minus_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - s;
}

inline void minus_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] - v2[i];
}

// v = v * s
void times_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
	for (Long i = 0; i < N; ++i) {
		v[i] = v1[i] * s;
	}
}

// v = v * v
inline void times_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] * v2[i];
}

// v = v / s
inline void divide_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
	T2 sInv{ 1 / s };
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] * sInv;
}

// v = s / v
inline void divide_vsv(Doub *v, Doub_I s, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = s / v1[i];
}

// v = v / v
inline void divide_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i] / v2[i];
}

// v = real(v)

inline void real_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = real(v1[i]); 
}

// v = imag(v)

inline void imag_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = imag(v1[i]);
}

inline void abs_vv(Doub *v, const Comp *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = abs(v1[i]);
}

inline void to_comp_vv(Comp *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = v1[i];
}

inline void dot_svv(Doub_O s, const Doub *v1, const Doub *v2, Long_I N)
{
	s = 0.;
	for (Long i = 0; i < N; ++i) {
		s += v1[i] * v2[i];
	}
}

inline void dot_svv(Comp_O s, const Comp *v1, const Comp *v2, Long_I N)
{
	s = 0.;
	for (Long i = 0; i < N; ++i) {
		s += conj(v1[i]) * v2[i];
	}
}

inline void dot_svv(Comp_O s, const Comp *v1, const Doub *v2, Long_I N)
{
	s = 0.;
	for (Long i = 0; i < N; ++i) {
		s += conj(v1[i]) * v2[i];
	}
}

inline void dot_svv(Comp_O s, const Doub *v1, const Comp *v2, Long_I N)
{
	s = 0.;
	for (Long i = 0; i < N; ++i) {
		s += v1[i] * v2[i];
	}
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

inline void invSqrt_vv(Float *v, const Float *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = 1.f / sqrt(v1[i]);
}

inline void invSqrt_vv(Doub *v, const Doub *v1, Long_I N)
{
	for (Long i = 0; i < N; ++i)
		v[i] = 1. / sqrt(v1[i]);
}

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

} // nemaspace slisc
