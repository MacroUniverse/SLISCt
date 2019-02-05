// meta-programming utilities

#pragma once
#include <type_traits>
#include "typedef.h"

namespace slisc {

// type traits

// declarations
template <class T, class U> struct promo_type;

// true_type() or true_type::value is true
// false_type() or false_type::value is false
using std::true_type; using std::false_type;

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
template <class T> struct is_complex : false_type {};
template <class T> struct is_complex<std::complex<T>> : true_type {};

// is_slisc<>
template <class T1, class T2 = void, class Enable = void> struct is_slisc : false_type {};

// single type param
// for T = FixVec<>, FixCmat<>, Vector<>, Matrix<>, Cmat<>, Diag<>, MatCoo<>, MatCooH<>
// is_slisc<T>() returns true, is_slisc<T>::type = value type of container
// otherwise is_slisc<T>() returns false, is_slisc<T>::type is undefined
template <class T> struct is_slisc<Vector<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Matrix<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Cmat<T>> : true_type { typedef T type; };
template <class T, Long Nr, Long Nc> struct is_slisc<FixVec<T, Nr, Nc>> : true_type { typedef T type; };
template <class T, Long Nr, Long Nc> struct is_slisc<FixCmat<T, Nr, Nc>> : true_type { typedef T type; };
template <class T> struct is_slisc<Mat3d<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Diag<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<MatCoo<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<MatCooH<T>> : true_type { typedef T type; };

// double type param
// equivalent to is_slisc<T1>() && is_slisc<T2>()
//template <class T1, class T2>
//struct is_slisc<T1, T2, typename std::conjunction<is_slisc<T1>, is_slisc<T2>>::type> : true_type
//{ typedef typename promo_type<T1, T2>::type type; };

// === type mapping ===

// map each scalr type to a unique number
// 0-19: integral types
// 20-39: floating point types
// 40-59: floating point complex types
// a smaller type in each level can be losslessly converted to a larger type in that level
// floating point type number plus 20 is the corresponding complex type number
template <class T>
constexpr Int type_num()
{
	if (is_integral<T>()) {
		if (is_same<T, Bool>()) return 0;
		if (is_same<T, Char>()) return 1;
		if (is_same<T, Uchar>()) return 2;
		if (is_same<T, Int>()) return 3;
		if (is_same<T, Uint>()) return 4;
		if (is_same<T, Llong>()) return 5;
		if (is_same<T, Ullong>()) return 6;
	}
	else if (is_floating_point<T>()) {
		if (is_same<T, float>()) return 20;
		if (is_same<T, Doub>()) return 21;
		if (is_same<T, Ldoub>()) return 22;
	}
	else if (is_complex<T>()) {
		if (is_same<T, std::complex<float>>()) return 40;
		if (is_same<T, Comp>()) return 41;
		if (is_same<T, std::complex<Ldoub>>()) return 42;
	}
	return -1;
}

// implementation of promo_type<T,U>
template <class T, class U>
constexpr auto promo_type0()
{
	static_assert(type_num<T>() != -1 && type_num<U>() != -1, "unhandled type!");

	constexpr Int Tnum = type_num<T>(), Unum = type_num<U>();
	if constexpr (Tnum < 40 && Unum < 40 || Tnum >= 40 && Unum >= 40) {
		// none are complex or both are complex, choose the larger type
		if constexpr (Tnum > Unum)
			return T();
		else
			return U();
	}
	else {
		// one complex one not
		if constexpr (Tnum < Unum) {
			// U is complex
			if constexpr (Unum - Tnum >= 20)
				return U(); // U is compatible
			else
				return std::complex<T>(); // need a larger complex
		}
		else {
			// T is complex
			if constexpr (Tnum - Unum >= 20)
				return T(); // T is compatible
			else
				return std::complex<U>(); // need a larger complex
		}
	}
}

// promo_type<T,U>::type is the smallest type that T, U can losslessly converted to
// (including from integer to floating point conversion)
// will be used as the return type of operator+-*/(T t, U u), etc.
// e.g. Bool + Doub = Doub; Int + Comp = Comp; Doub + complex<float> = Comp;
// see promo_type0() for detail
template <class T, class U> struct promo_type
{ using type = decltype(promo_type0<T, U>()); };

// for T = bool, character types
// to_num_t<T>::type is Int
// otherwise is T
template <typename T> struct num_type { typedef T type; };
template<> struct num_type<Bool> { typedef Int type; };
template<> struct num_type<Char> { typedef Int type; };
template<> struct num_type<Uchar> { typedef Int type; };

// for T = std::complex<T1>
// rm_comp<T>::type is T1
// otherwise is T
template <class T> struct rm_complex { typedef T type; };
template <class T> struct rm_complex<std::complex<T>> { typedef T type; };

}
