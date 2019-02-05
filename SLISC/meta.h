// meta-programming utilities

#pragma once
#include <type_traits>
#include "typedef.h"

namespace slisc {

// type traits

// for U = T
// is_same<U,T>() or is_same<U,T>::value returns true
// otherwise returns false
using std::is_same;

// for T = c++ primitive types
// is_fundamental<T>() or is_fundamental<T>::value returns true
// otherwise returns false
using std::is_fundamental;

// for T = float, double or long double,
// is_floating_point<T>() or is_floating_point<T>::value returns true
// otherwise returns false
using std::is_floating_point;

// for T = bool, character types, integer types
// is_integral<T>() or is_integral<T>::value returns true
// otherwise returns false
using std::is_integral;

// for T = std::complex<T1>
// is_complex<T>() or is_complex<T>::value returns true
// otherwise returns false
template <class T> struct is_complex : std::false_type {};
template <class T> struct is_complex<std::complex<T>> : std::true_type {};

// === type mapping ===

// map each type to a unique number
// a smaller type can be losslessly converted to a larger type
// or can be converted from an integral type to floating point type
// or from a primitive type to a floating point complex type
template <class T>
constexpr Int type_num()
{
	if (is_integral<T>()) {
		if (is_same<T, Bool>()) return 1;
		if (is_same<T, Char>()) return 2;
		if (is_same<T, Uchar>()) return 3;
		if (is_same<T, Int>()) return 4;
		if (is_same<T, Uint>()) return 5;
		if (is_same<T, Llong>()) return 6;
		if (is_same<T, Ullong>()) return 7;
	}
	else if (is_floating_point<T>()) {
		if (is_same<T, float>()) return 11;
		if (is_same<T, Doub>()) return 12;
		if (is_same<T, Ldoub>()) return 13;
	}
	else if (is_complex<T>()) {
		if (is_same<T, std::complex<float>>()) return 21;
		if (is_same<T, Comp>()) return 22;
		if (is_same<T, std::complex<Ldoub>>()) return 23;
	}
	return -1;
}

// return the default value of the larger type
template <class T, class U>
constexpr auto promo_type0()
{
	static_assert(type_num<T> != -1 && type_num<U> != -1, "unhandled type!");
	if constexpr (type_num<T>() > type_num<U>())
		return T();
	else
		return U();
}

// promo_type<T, U>::type is the larger type of T and U
template <class T, class U> struct promo_type
{
	using type = decltype(promo_type0<T, U>());
};

// for T = bool, character types
// to_num_t<T>::type is Int
// otherwise is T
template <typename T> struct to_num_t { typedef T type; };
template<> struct to_num_t<Bool> { typedef Int type; };
template<> struct to_num_t<Char> { typedef Int type; };
template<> struct to_num_t<Uchar> { typedef Int type; };

// for T = std::complex<T1>
// rm_comp<T>::type is T1
// otherwise is T
template <class T> struct rm_comp { typedef T type; };
template <class T> struct rm_comp<std::complex<T>> { typedef T type; };

}
