// scalar utilities
#pragma once
#include "meta.h"

namespace slisc {

// check if two scalars have the same types and values
// const-ness and reference are ignored
template <class T1, class T2>
Bool is_equiv(const T1 &s1, const T2 &s2)
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
	return (const typename to_num_t<T>::type)s;
}

// basic functions
using std::abs; using std::sqrt;
using std::sin; using std::cos; using std::tan;

template <class T>
inline auto cot(const T &x)
{ return 1 / tan(x); }

// operator+,-,*,/ between std::complex<intrinsic types> and intrinsic types

template <class T, class Tc>
auto operator+(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_fundamental<T>() && is_fundamental<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z) + x, imag(z));
}

template <class T, class Tc>
auto operator+(const T &x, const std::complex<Tc> &z)
{
	return z + x;
}

template <class T, class Tc>
auto operator-(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_fundamental<T>() && is_fundamental<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z) - x, imag(z));
}

template <class T, class Tc>
auto operator-(const T &x, const std::complex<Tc> &z)
{
	static_assert(is_fundamental<T>() && is_fundamental<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(x - real(z), -imag(z));
}

template <class T, class Tc>
auto operator*(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_fundamental<T>() && is_fundamental<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z)*x, imag(z)*x);
}

template <class T, class Tc>
auto operator*(const T &x, const std::complex<Tc> &z)
{
	return z * x;
}

template <class T, class Tc>
auto operator/(const std::complex<Tc> &z, const T &x)
{
	static_assert(is_fundamental<T>() && is_fundamental<Tc>(), "type error!");
	return std::complex<typename promo_type<T, Tc>::type>(real(z)/x, imag(z)/x);
}

template <class T, class Tc>
auto operator/(const T &x, const std::complex<Tc> &z)
{ return (1./z) * x; }

// operator+,-,*,/ between two different std::complex<> types

template <class Tx, class Ty>
auto operator+(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_fundamental<Tx>() && is_fundamental<Ty>(), "type error!");
	return std::complex<typename promo_type<Tx, Ty>::type>
		(real(x) + real(y), imag(x) + imag(y));
}

template <class Tx, class Ty>
auto operator-(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_fundamental<Tx>() && is_fundamental<Ty>(), "type error!");
	return std::complex<typename promo_type<Tx, Ty>::type>
		(real(x) - real(y), imag(x) - imag(y));
}

// TODO: slow
template <class Tx, class Ty>
auto operator*(const std::complex<Tx> &x, const std::complex<Ty> &y)
{
	static_assert(is_fundamental<Tx>() && is_fundamental<Ty>(), "type error!");
	return std::complex<typename promo_type<Tx, Ty>::type> 
		(real(x)*real(y) - imag(x)*imag(y), real(x)*imag(y) + imag(x)*real(y));
}

// TODO: slow
template <class Tx, class Ty>
auto operator/(const std::complex<Tx> &x, const std::complex<Ty> &y)
{ return x * (1 / y); }

}
