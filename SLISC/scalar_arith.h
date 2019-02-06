// scalar utilities
#pragma once
#include <cmath>
#include <algorithm>
#include "meta.h"

namespace slisc {

// constant value of a type

// CONST<T, s>() will generate a non-complex constexpr value even if T is complex
template <class T, Int val>
constexpr const auto CONST()
{
	typedef typename rm_complex<T>::type Real;
	constexpr Real ret(val);
	return ret;
}

// CONST<T, real, imag>(), where T must be std::complex<>, will generate a constexpr complex value
template <class T, Int real, Int imag, class Enable = typename std::enable_if<is_complex<T>()>::type>
constexpr const T CONSTC()
{
	constexpr T ret(real, imag);
	return ret;
}

// these are slightly different from Numerical Recipes
template<class T1, class T2>
constexpr const auto &MIN(const T1 &a, const T2 &b)
{ return a < b ? a : b; }

template<class T1, class T2>
constexpr const auto &MAX(const T1 &a, const T2 &b)
{ return a < b ? b : a; }

template<class T>
constexpr const T SQR(const T &a) { return a * a; }

template<class T>
constexpr const T SIGN(const T &s)
{
	constexpr T one = T(1), neg_one = T(-1), zero = T();
	return 0 < s ? one : (s < 0 ? neg_one : zero);
}

template<class T1, class T2>
constexpr const T1 SIGN(const T1 &a, const T2 &b)
{ return 0 <= b ? (0 <= a ? a : -a) : (0 <= a ? -a : a); }

template<class T>
inline void SWAP(T &a, T &b)
{ error("use std::swap() instead!"); }

// check if two scalars have the same types and values
// const-ness and reference are ignored
template <class T1, class T2>
constexpr Bool is_equiv(const T1 &s1, const T2 &s2)
{
	if constexpr (is_same<T1, T2>()) {
		return s1 == s2;
	}
	return false;
}

// convert false to 0, true to 1, character to Int
template <class T>
auto to_num(const T &s)
{
	return (const typename num_type<T>::type)s;
}

// basic functions
using std::swap; using std::abs; using std::real; using std::imag;
using std::sqrt; using std::sin; using std::cos; using std::tan;
using std::exp; using std::log; using std::log10;
using std::expm1; using std::log1p; using std::hypot;
using std::sinh; using std::cosh; using std::tanh;

template <class T>
inline auto cot(const T &x)
{
	typedef typename rm_complex<T>::type Tr;
	constexpr Tr one = Tr(1);
	return one / tan(x);
}

template <class T>
inline T sinc(const T &x)
{
	typedef typename rm_complex<T>::type Tr;
	constexpr Tr one = Tr(1);
	return x == 0 ? one : sin(x) / x;
}

// check if an integer is odd
inline Bool isodd(Int_I n) { return n & 1; }
inline Bool isodd(Long_I n) { return n & 1; }

// return true if n is power of 2 or 0
inline Bool ispow2(Int_I n) { return (n&(n-1)) == 0; }
inline Bool ispow2(Long_I n) { return (n&(n-1)) == 0; }

// return the positive modulus (use "%" when i >= 0)
inline Int mod(Int_I i, Int_I n) { return (i % n + n) % n; }
inline Long mod(Long_I i, Long_I n) { return (i % n + n) % n; }

// matrix double index to single index conversion

inline Long csub2ind(Long_I Nr, Long_I i, Long_I j)
{ return i + Nr*j; } // column major

inline Long rsub2ind(Long_I Nc, Long_I i, Long_I j)
{ return Nc*i + j; } // row major


// operator+,-,*,/ between floating point std::complex<> and all arithmetic types
// return promo_type<>

template <class T, class Tc>
const auto operator+(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_arithmetic<T>() && is_floating_point<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z) + x, imag(z));
}

template <class T, class Tc>
const auto operator+(const T &x, const std::complex<Tc> &z)
{ return z + x; }

template <class T, class Tc>
const auto operator-(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_arithmetic<T>() && is_floating_point<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z) - x, imag(z));
}

template <class T, class Tc>
const auto operator-(const T &x, const std::complex<Tc> &z)
{
	static_assert(is_arithmetic<T>() && is_floating_point<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(x - real(z), -imag(z));
}

template <class T, class Tc>
const auto operator*(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_arithmetic<T>() && is_floating_point<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z)*x, imag(z)*x);
}

template <class T, class Tc>
const auto operator*(const T &x, const std::complex<Tc> &z)
{ return z * x; }

template <class T, class Tc>
const auto operator/(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_arithmetic<T>() && is_floating_point<Tc>(), "type error!");
	typedef typename promo_type<T, Tc>::type Tp;
	Tp inv_x; // Tp is floating point
	inv_x = Tp(1) / x; // T is integral
	return z * inv_x;
}

template <class T, class Tc>
const auto operator/(const T &x, const std::complex<Tc> &z)
{ return x * (Tc(1)/z); }

// operator+,-,*,/ between two different std::complex<> types
// return promo type

template <class Tx, class Ty>
const auto operator+(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_floating_point<Tx>() && is_floating_point<Ty>(), "type error!");
	return std::complex<typename promo_type<Tx, Ty>::type>
		(real(x) + real(y), imag(x) + imag(y));
}

template <class Tx, class Ty>
const auto operator-(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_floating_point<Tx>() && is_floating_point<Ty>(), "type error!");
	return std::complex<typename promo_type<Tx, Ty>::type>
		(real(x) - real(y), imag(x) - imag(y));
}

template <class Tx, class Ty>
const auto operator*(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_floating_point<Tx>() && is_floating_point<Ty>(), "type error!");
	typedef std::complex<typename promo_type<Tx, Ty>::type> T;
	return T(x) * T(y);
}

template <class Tx, class Ty>
const auto operator/(const std::complex<Tx> &x, const std::complex<Ty> &y)
{ return x * (Ty(1) / y); }

}