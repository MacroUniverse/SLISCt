// meta-programming utilities

#pragma once
#include <type_traits>
#include "global.h"

// must use SLISC_IF(()) 
#define SLISC_IF(cond) typename std::enable_if<cond>::type* = 0

namespace slisc {

// type traits

// declarations
template <class T, class U> struct promo_type;
template <class T> struct rm_complex;
template <class T> constexpr Int type_num();

// true_type::value is true, false_type::value is false
// true_type() and false_type() can be implicitly converted to true and false
using std::true_type; using std::false_type;

// integral_constant<Bool, true> is true_type
// integral_constant<Bool, false> is false_type
using std::integral_constant;

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

using std::is_signed;

// for T = floating point std::complex<>
// is_fp_complex<T>() or is_fp_complex<T>::value returns true
// otherwise returns false
template<class T> struct is_fp_comp : false_type {};

template<class T>
struct is_fp_comp<std::complex<T>> :
	integral_constant<bool, is_floating_point<T>::value> {};

// map each scalr type to a unique number
// 0-19: integral types
// 20-39: floating point types
// 40-59: floating point complex types
// a smaller type in each level can be losslessly converted to a larger type in that level
// floating point type number plus 20 is the corresponding complex type number
template <class T>
constexpr Int type_num()
{
	if constexpr (is_integral<T>()) {
		if (is_same<T, Bool>()) return 0;
		if (is_same<T, Char>()) return 1;
		if (is_same<T, Int>()) return 2;
		if (is_same<T, Llong>()) return 3;
	}
	else if constexpr (is_floating_point<T>()) {
		if (is_same<T, float>()) return 20;
		if (is_same<T, Doub>()) return 21;
		if (is_same<T, Ldoub>()) return 22;
	}
	else if constexpr (is_fp_comp<T>()) {
		if (is_same<T, std::complex<float>>()) return 40;
		if (is_same<T, Comp>()) return 41;
		if (is_same<T, std::complex<Ldoub>>()) return 42;
	}
	return -1;
}

// check if is listed scalar
template<class T>
struct is_scalar : integral_constant<bool,
	type_num<T>() >= 0> {};

template<class T>
struct is_real : integral_constant<bool,
	type_num<T>() >= 0 && !is_fp_comp<T>()> {};

// check if is specific container type
template <class T> struct is_Vector : false_type {};
template <class T> struct is_Vector<Vector<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_Matrix : false_type {};
template <class T> struct is_Matrix<Matrix<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_Cmat : false_type {};
template <class T> struct is_Cmat<Cmat<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_FixVec : false_type {};
template <class T, Long N> struct is_FixVec<FixVec<T, N>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_FixCmat : false_type {};
template <class T, Long Nr, Long Nc> struct is_FixCmat<FixCmat<T, Nr, Nc>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_Mat3d : false_type {};
template <class T> struct is_Mat3d<Mat3d<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_Diag : false_type {};
template <class T> struct is_Diag<Diag<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_MatCoo : false_type {};
template <class T> struct is_MatCoo<MatCoo<T>> : integral_constant<Bool, is_scalar<T>::value> {};

template <class T> struct is_MatCooH : false_type {};
template <class T> struct is_MatCooH<MatCooH<T>> : integral_constant<Bool, is_scalar<T>::value> {};

// check if is fixed-size container
template <class T> struct is_fixed : integral_constant<Bool,
	is_FixVec<T>() || is_FixCmat<T>()>{};

// check if is dense container (including fixed-size)
template <class T> struct is_dense : integral_constant<Bool,
	is_Vector<T>() || is_Matrix<T>() || is_Cmat<T>() || is_fixed<T>()
	|| is_Mat3d<T>()>{};

// check if is sparse vector/matrix
template <class T> struct is_sparse : integral_constant<Bool,
	is_Diag<T>() || is_MatCoo<T>() || is_MatCooH<T>()> {};

// check if is any slisc container
template <class T> struct is_contain : integral_constant<Bool,
	is_dense<T>() || is_sparse<T>()> {};

template <class T>
constexpr Int contain_num()
{
	if constexpr (is_Vector<T>()) return 0;
	else if constexpr (is_Matrix<T>()) return 1;
	else if constexpr (is_Cmat<T>()) return 2;
	else if constexpr (is_Mat3d<T>()) return 3;

	else if constexpr (is_FixVec<T>()) return 20;
	else if constexpr (is_FixCmat<T>()) return 22;

	else if constexpr (is_Diag<T>()) return 31;
	else if constexpr (is_MatCoo<T>()) return 32;
	else if constexpr (is_MatCooH<T>()) return 33;

	return -1;
}

// check if two containers are the same (value_type can be different)
template <class T1, class T2> struct is_same_contain : integral_constant<Bool,
	contain_num<T1>() >= 0 && contain_num<T1>() == contain_num<T2>()> {};

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
template <class T, Int val, SLISC_IF((is_arithmetic<typename rm_complex<T>::type>::value))>
struct Aconst
{
	static constexpr typename rm_complex<T>::type value = val;
};

// Cconst<T, real, imag>::value, where T must be std::complex<>, is a constexpr variable of complex<>
// will not instanciate if T is illegal
template <class T, Int real, Int imag = 0, SLISC_IF((is_fp_comp<T>::value))>
struct Cconst
{
	static constexpr T value = T(real, imag);
};

// === type mapping ===

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
