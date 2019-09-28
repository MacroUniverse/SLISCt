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

inline void plus_equals_vs(Int *v, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] += s;
}

inline void plus_equals_vs(Long *v, Long_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] += s;
}

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

inline void plus_equals_vs(Comp *v, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] += s;
}

// v += v

inline void plus_equals_vv(Int *v, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] += v1[i];
}

inline void plus_equals_vv(Long *v, const Long *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] += v1[i];
}

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

// v -= s

inline void minus_equals_vs(Char *v, Char_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Int *v, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Long *v, Long_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Float *v, Float_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Doub *v, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Comp *v, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

inline void minus_equals_vs(Comp *v, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= s;
}

// v -= v

inline void minus_equals_vv(Char *v, const Char *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

inline void minus_equals_vv(Int *v, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

inline void minus_equals_vv(Float *v, const Float *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

inline void minus_equals_vv(Doub *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

inline void minus_equals_vv(Comp *v, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

inline void minus_equals_vv(Comp *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] -= v1[i];
}

// v *= s

inline void times_equals_vs(Int *v, Int_I s, Long N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= s;
}

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

// v *= v

inline void times_equals_vv(Char *v, const Char *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

inline void times_equals_vv(Int *v, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

inline void times_equals_vv(Float *v, const Float *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

inline void times_equals_vv(Doub *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

inline void times_equals_vv(Comp *v, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

inline void times_equals_vv(Comp *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= v1[i];
}

// v /= s

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

inline void divide_equals_vs(Float *v, const Float &s, Long_I N)
{ times_equals_vs(v, 1.f/s, N); }

inline void divide_equals_vs(Doub *v, const Doub &s, Long_I N)
{ times_equals_vs(v, 1./s, N); }

inline void divide_equals_vs(Comp *v, const Doub &s, Long_I N)
{ times_equals_vs(v, 1./s, N); }

// v /= v

inline void divide_equals_vv(Int *v, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] /= v1[i];
}

inline void divide_equals_vv(Float *v, const Float *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] /= v1[i];
}

inline void divide_equals_vv(Doub *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] /= v1[i];
}

inline void divide_equals_vv(Comp *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] /= v1[i];
}

inline void divide_equals_vv(Comp *v, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] /= v1[i];
}

// v = mod(v, s)

inline void mod_vvs(Char *v, const Char *v1, Char_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = mod(v1[i], s);
}

inline void mod_vvs(Int *v, const Int *v1, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = mod(v1[i], s);
}

inline void mod_vvs(Long *v, const Long *v1, Long_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = mod(v1[i], s);
}

// v %= s

inline void rem_vs(Int *v, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] %= s;
}

inline void rem_vs(Long *v, Long_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] %= s;
}

// v = v % s

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

// v = v + s

inline void plus_vvs(Char *v, const Char *v1, Char_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Int *v, const Int *v1, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Float *v, const Float *v1, Float_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Comp *v, const Doub *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Comp *v, const Comp *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

inline void plus_vvs(Comp *v, const Comp *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + s;
}

// v = v + v

inline void plus_vvv(Int *v, const Int *v1, const Int *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + v2[i];
}

inline void plus_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + v2[i];
}

inline void plus_vvv(Comp *v, const Comp *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] + v2[i];
}

// -v
inline void minus_v(Int *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= -1;
}

inline void minus_v(Float *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= -1;
}

inline void minus_v(Doub *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= -1;
}

inline void minus_v(Comp *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] *= -1;
}


// v = -v

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

inline void minus_vv(Comp *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = -v1[i];
}

inline void minus_vv(Comp *v, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = -v1[i];
}

// v = s - v

inline void minus_vsv(Int *v, Int_I s, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s - v1[i];
}

inline void minus_vsv(Doub *v, Doub_I s, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s - v1[i];
}

inline void minus_vsv(Comp *v, Doub_I s, const Comp *v1, Long_I N)
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

inline void minus_vvs(Float *v, const Float *v1, Float_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - s;
}

inline void minus_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - s;
}

inline void minus_vvs(Comp *v, const Doub *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - s;
}

inline void minus_vvs(Comp *v, const Comp *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - s;
}

inline void minus_vvs(Comp *v, const Comp *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - s;
}

// v = v - v

inline void minus_vvv(Int *v, const Int *v1, const Int *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - v2[i];
}

inline void minus_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - v2[i];
}

inline void minus_vvv(Comp *v, const Comp *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - v2[i];
}

inline void minus_vvv(Comp *v, const Doub *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - v2[i];
}

inline void minus_vvv(Comp *v, const Comp *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] - v2[i];
}

// v = v * s

void times_vvs(Int *v, const Int *v1, Int_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

void times_vvs(Float *v, const Float *v1, Float_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

void times_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

void times_vvs(Comp *v, const Doub *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

void times_vvs(Comp *v, const Comp *v1, Doub_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

void times_vvs(Comp *v, const Comp *v1, Comp_I s, Long_I N)
{
    for (Long i = 0; i < N; ++i) {
        v[i] = v1[i] * s;
    }
}

// v = v * v

inline void times_vvv(Int *v, const Int *v1, const Int *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

inline void times_vvv(Float *v, const Float *v1, const Float *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

inline void times_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

inline void times_vvv(Comp *v, const Doub *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

inline void times_vvv(Comp *v, const Comp *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

inline void times_vvv(Comp *v, const Comp *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] * v2[i];
}

// v = v / s

inline void divide_vvs(Float *v, const Float *v1, Float_I s, Long_I N)
{ times_vvs(v, v1, 1.f / s, N); }

inline void divide_vvs(Doub *v, const Doub *v1, Doub_I s, Long_I N)
{ times_vvs(v, v1, 1. / s, N); }

inline void divide_vvs(Comp *v, const Comp *v1, Doub_I s, Long_I N)
{ times_vvs(v, v1, 1. / s, N); }

inline void divide_vvs(Comp *v, const Doub *v1, Comp_I s, Long_I N)
{ times_vvs(v, v1, 1. / s, N); }

inline void divide_vvs(Comp *v, const Comp *v1, Comp_I s, Long_I N)
{ times_vvs(v, v1, 1. / s, N); }

// v = s / v

inline void divide_vsv(Float *v, Float_I s, const Float *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s / v1[i];
}

inline void divide_vsv(Doub *v, Doub_I s, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s / v1[i];
}

inline void divide_vsv(Comp *v, Doub_I s, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s / v1[i];
}

inline void divide_vsv(Comp *v, Comp_I s, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s / v1[i];
}

inline void divide_vsv(Comp *v, Comp_I s, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = s / v1[i];
}

// v = v / v

inline void divide_vvv(Int *v, const Int *v1, const Int *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] / v2[i];
}

inline void divide_vvv(Doub *v, const Doub *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] / v2[i];
}

inline void divide_vvv(Comp *v, const Doub *v1, const Comp *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] / v2[i];
}

inline void divide_vvv(Comp *v, const Comp *v1, const Doub *v2, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i] / v2[i];
}

inline void divide_vvv(Comp *v, const Comp *v1, const Comp *v2, Long_I N)
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

// real(v)

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

inline void abs_v(Int *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v[i]);
}

inline void abs_v(Doub *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v[i]);
}

inline void abs_v(Comp *v, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v[i]);
}

// v = abs(v)

inline void abs_vv(Int *v, const Int *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v1[i]);
}

inline void abs_vv(Doub *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v1[i]);
}

inline void abs_vv(Doub *v, const Comp *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = abs(v1[i]);
}

// v = comp(v)

inline void to_comp_vv(Comp *v, const Float *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i];
}

inline void to_comp_vv(Comp *v, const Doub *v1, Long_I N)
{
    for (Long i = 0; i < N; ++i)
        v[i] = v1[i];
}

// s = sum(v)

inline Long sum_v(const Bool *v, Long_I N)
{
    Long s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

inline Int sum_v(const Int *v, Long_I N)
{
    Int s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

inline Long sum_v(const Long *v, Long_I N)
{
    Long s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

inline Float sum_v(const Float *v, Long_I N)
{
    Float s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

inline Doub sum_v(const Doub *v, Long_I N)
{
    Doub s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

inline Comp sum_v(const Comp *v, Long_I N)
{
    Comp s = v[0];
    for (Long i = 1; i < N; ++i)
        s += v[i];
    return s;
}

// s = max(v)

inline Int max_v(const Int *v, Long_I N)
{
    Int s = v[0], val;
    for (Long i = 1; i < N; ++i) {
        if (s < v[i])
            s = v[i];
    }
    return s;
}

inline Float max_v(const Float *v, Long_I N)
{
    Float s = v[0], val;
    for (Long i = 1; i < N; ++i) {
        if (s < v[i])
            s = v[i];
    }
    return s;
}

inline Doub max_v(const Doub *v, Long_I N)
{
    Doub s = v[0], val;
    for (Long i = 1; i < N; ++i) {
        if (s < v[i])
            s = v[i];
    }
    return s;
}

// s = max_abs(v)

inline Doub max_abs(const Doub *v, Long_I N)
{
    Doub s = abs(v[0]), val;
    for (Long i = 1; i < N; ++i) {
        val = abs(v[i]);
        if (s < val)
            s = val;
    }
    return s;
}

inline Doub max_abs_v(const Comp *v, Long_I N)
{
    Doub s = abs(v[0]), val;
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
