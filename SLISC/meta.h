// meta-programming utilities

#pragma once
#include <type_traits>
#include "typedef.h"

namespace slisc {

// type traits

// declarations
template <class T, class U> struct promo_type;
template <class T> struct rm_complex;

// true_type() or true_type::value is true
// false_type() or false_type::value is false
using std::true_type; using std::false_type;

// for U = T
// is_same<U,T>() or is_same<U,T>::value returns true
// otherwise returns false
using std::is_same;

// for T = bool, character types, integer types
// is_integral<T>() or is_integral<T>::value returns true
// otherwise returns false
using std::is_integral;

// for T = float, double, long double,
// is_floating_point<T>() or is_floating_point<T>::value returns true
// otherwise returns false
using std::is_floating_point;

// for T = integral type, floating point type
// is_arithmetic<T>() or is_arithmetic<T>::value returns true
// otherwise returns false
using std::is_arithmetic;

// for T = arithmetic type, void, null pointer
// is_fundamental<T>() or is_fundamental<T>::value returns true
// otherwise returns false
using std::is_fundamental;

// for T = arithmetic std::complex<>
// is_complex<T>() or is_complex<T>::value returns true
// otherwise returns false
template<class T> struct is_complex : false_type {};

template<class T>
struct is_complex<std::complex<T>> :
	integral_constant<bool, is_arithmetic<T>::value> {};

// for T = arithmetic type, complex<> type
// is_scalar<T>() or is_scalar<T>::value returns true
// otherwise returns false
template< class T >
struct is_scalar : std::integral_constant<bool,
	is_arithmetic<T>() || is_complex<T>()
> {};

// is_slisc<>
template <class T> struct is_slisc : false_type {};

// single type param
// for T = FixVec<>, FixCmat<>, Vector<>, Matrix<>, Cmat<>, Diag<>, MatCoo<>, MatCooH<>
// is_slisc<T>() returns true, is_slisc<T>::type = value type of container
// otherwise is_slisc<T>() returns false, is_slisc<T>::type is undefined
template <class T> struct is_slisc<Vector<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Matrix<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Cmat<T>> : true_type { typedef T type; };
template <class T, Long N> struct is_slisc<FixVec<T, N>> : true_type { typedef T type; };
template <class T, Long Nr, Long Nc> struct is_slisc<FixCmat<T, Nr, Nc>> : true_type { typedef T type; };
template <class T> struct is_slisc<Mat3d<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<Diag<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<MatCoo<T>> : true_type { typedef T type; };
template <class T> struct is_slisc<MatCooH<T>> : true_type { typedef T type; };

// double type param
// if (is_slisc<T1>() && is_slisc<T2>())
template <class T1, class T2, class Enable = void> struct is_slisc2 : false_type {};

template <class T1, class T2>
struct is_slisc2<T1, T2,
	typename std::enable_if<
		std::conjunction<is_slisc<T1>, is_slisc<T2>>()
	>::type
> : true_type {
	typedef typename
		promo_type<typename T1::value_type, typename T2::value_type>::type
		type;
};

// === get static constexpr values ===

// Sconst<T, val>::value is a constexpr variable of arithmetic type or complex<>
template <class T, Int val>
struct Sconst
{
	static constexpr T value = val;
};

// Aconst<T, val>::value is a constexpr variable of arithmetic type
// when T is complex<value_type>, CONST<T, val>::value is value_type
// will not instanciate if T is illegal
template <class T, Int val, class Enable =
	typename enable_if<is_arithmetic<typename rm_complex<T>::type>()>::type>
struct Aconst
{
	static constexpr typename rm_complex<T>::type value = val;
};

// Cconst<T, real, imag>::value, where T must be std::complex<>, is a constexpr variable of complex<>
// will not instanciate if T is illegal
template <class T, Int real, Int imag = 0, class Enable = typename std::enable_if<is_complex<T>()>::type>
struct Cconst
{
	static constexpr T value = T(real, imag);
};

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
		if (is_same<T, Int>()) return 2;
		if (is_same<T, Llong>()) return 3;
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
	error("unsupported type!");
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
				return U(); // U is lossless
			else
				return std::complex<T>(); // need a larger complex
		}
		else {
			// T is complex
			if constexpr (Tnum - Unum >= 20)
				return T(); // T is lossless
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
