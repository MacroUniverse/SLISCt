// scalar utilities
#pragma once
#include "complex_arith.h"
#include "imag.h"

namespace slisc {

// these are slightly different from Numerical Recipes

template <class T, SLS_IF(is_scalar<T>())>
inline Bool ISNAN(const T &s)
{
    return s != s;
}

template<class T>
constexpr const T &MIN(const T &a, const T &b)
{ return a < b ? a : b; }

template<class T>
constexpr const T &MAX(const T &a, const T &b)
{
    return a < b ? b : a;
}

template<class T1, class T2, SLS_IF(!is_same<T1, T2>())>
constexpr auto MAX(const T1 &a, const T2 &b)
{
    return a < b ? b : a;
}

template<class T>
constexpr const T SQR(const T &a) { return a * a; }

template<class T, SLS_IF(is_fpt<T>())>
constexpr const T ABS2(const T &a)
{
    return a * a;
}

template<class T, SLS_IF(is_comp<T>())>
constexpr const rm_comp<T> ABS2(const T &a)
{
    return SQR(real(a)) + SQR(imag(a));
}

constexpr Float SIGN(Float_I s)
{ return s > 0.f ? 1.f : (s < 0.f ? -1.f : 0.f); }

constexpr Doub SIGN(Doub_I s)
{ return s > 0. ? 1. : (s < 0. ? -1. : 0.); }

template<class T>
inline T SIGN(const T &a, const T &b)
{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

// get std::vector size in Long instead of size_t
template <class T, SLS_IF(is_vector<T>() || is_basic_str<T>())>
inline Long Size(const T &v)
{
    return (Long)v.size();
}

// for SWAP, use std::swap instead

template <class T, SLS_IF(is_fpt<T>())>
inline T INV(const T &x)
{ return 1 / x; }

inline Fcomp INV(Fcomp_I x)
{ return 1.f/x; }

inline Comp INV(Comp_I x)
{ return 1./x; }

inline Lcomp INV(Lcomp_I x)
{ return 1.l/x; }

inline Comp CONJ(Comp_I x)
{ return conj(x); }

inline Doub CONJ(Doub_I x)
{ return x; }

// check if `elm` is one of `vec[i]`
// `vec.size()` must be defined
// if you need i, see `search()` in `search.h`
template <class T1, class T2>
inline Bool is_in(const T1 &elm, const T2 &vec)
{
    for (Long i = 0; i < Size(vec); ++i)
        if (elm == vec[i])
            return true;
    return false;
}

// check if two scalars have the same types and values
// const-ness and reference are ignored
template <class T1, class T2>
constexpr Bool is_equiv(const T1 &s1, const T2 &s2)
{
    if (is_same<T1, T2>()) {
        return s1 == s2;
    }
    return false;
}

// convert bool and character to Int, others unchanged
inline Int to_num(Bool_I x) { return (Int)x; }
inline Int to_num(Char_I x) { return (Int)x; }
inline Int to_num(Uchar_I x) { return (Int)x; }
template <class T>
inline const T to_num(const T &x) { return x; }

// integer functions
// check if an integer is odd
template <class Tint, SLS_IF(is_integral<Tint>())>
inline Bool isodd(const Tint &n)
{
    return n & 1;
}

// return true if n is power of 2 or 0
template <class Tint, SLS_IF(is_integral<Tint>())>
inline Bool ispow2(const Tint &n)
{
    return (n&(n-1)) == 0;
}

// return the positive modulus (use "%" when i >= 0)
template <class T1, class T2, SLS_IF(
    is_integral<T1>() && is_scalar<T1>() && is_promo<T1, T2>())>
inline T1 mod(const T1 &i, const T2 &n)
{
    return (i % n + n) % n;
}

// return the positive modulus (s = n * d + return)
template <class T1, class T2, SLS_IF(
    is_fpt<T1>() && is_scalar<T1>() && is_promo<T1, T2>())>
    inline T1 mod(Long_O n, const T1 &s, const T2 &d)
{
    n = floor(s/d);
    return s - n * d;
}

// matrix double index to single index conversion

inline Long csub2ind(Long_I Nr, Long_I i, Long_I j)
{ return i + Nr*j; } // column major

inline Long rsub2ind(Long_I Nc, Long_I i, Long_I j)
{ return Nc*i + j; } // row major

// floating point functions

using std::min; using std::max; using std::swap;
using std::abs; using std::real; using std::imag;
using std::conj; using std::pow;
using std::sqrt; using std::sin; using std::cos; using std::tan;
using std::exp; using std::log; using std::log10;
using std::expm1; using std::log1p; using std::hypot;
using std::sinh; using std::cosh; using std::tanh;

inline Float cot(Float_I x) { return 1.f / tan(x); }
inline Doub cot(Doub_I x) { return 1. / tan(x); }

inline Float sinc(Float_I x) { return x == 0.f ? 1.f : sin(x) / x; }
inline Doub sinc(Doub_I x) { return x == 0. ? 1. : sin(x) / x; }

// factorial of a number

inline Doub factorial_imp(Doub_I n) {
    if (n == 0. || n == 1.)
        return 1.;
    else
        return n * factorial_imp(n - 1.);
}

inline Doub factorial(Long_I n) {
    if (n > 150)
        SLS_ERR("n too large!");
    return factorial_imp((Doub)n);
}

} // namespace slisc
