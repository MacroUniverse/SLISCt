// meta-programming utilities

#pragma once
#include <type_traits>
#include "global.h"

// using variable argument macro to allow parsing of ",".
// otherwise, "," will separate a single argument into multiple arguments
#define SLS_IF_HELPER(cond) typename std::enable_if<(bool)(cond)>::type* = 0
#define SLS_IF(...) SLS_IF_HELPER(__VA_ARGS__)

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
struct is_fp_comp<complex<T>> :
	integral_constant<bool, is_floating_point<T>::value> {};

template<class T>
struct is_Bool : integral_constant<Bool, is_same<T, Bool>()> {};

template<class T>
struct is_Char : integral_constant<Bool, is_same<T, Char>()> {};


template<class T>
struct is_Int : integral_constant<Bool, is_same<T, Int>()> {};

template<class T>
struct is_Llong : integral_constant<Bool, is_same<T, Llong>()> {};

template<class T>
struct is_Float : integral_constant<Bool, is_same<T, Float>()> {};

template<class T>
struct is_Doub : integral_constant<Bool, is_same<T, Doub>()> {};

template<class T>
struct is_Ldoub : integral_constant<Bool, is_same<T, Ldoub>()> {};

template<class T>
struct is_Fcomp : integral_constant<Bool, is_same<T, Fcomp>()> {};

template<class T>
struct is_Comp : integral_constant<Bool, is_same<T, Comp>()> {};

template<class T>
struct is_Lcomp : integral_constant<Bool, is_same<T, Lcomp>()> {};

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
		if (is_Bool<T>()) return 0;
		if (is_Char<T>()) return 1;
		if (is_Int<T>()) return 2;
		if (is_Llong<T>()) return 3;
	}
	else if (is_floating_point<T>()) {
		if (is_Float<T>()) return 20;
		if (is_Doub<T>()) return 21;
		if (is_Ldoub<T>()) return 22;
	}
	else if (is_fp_comp<T>()) {
		if (is_Fcomp<T>()) return 40;
		if (is_Comp<T>()) return 41;
		if (is_Lcomp<T>()) return 42;
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

template <class T> struct is_dense_mat : integral_constant<Bool,
	is_Matrix<T>() || is_Cmat<T>() || is_FixCmat<T>()
	|| is_Mat3d<T>()> {};

// check if is dense container (including fixed-size)
template <class T> struct is_dense : integral_constant<Bool,
	is_Vector<T>() || is_FixVec<T>() || is_dense_mat<T>()>{};

template <class T> struct is_real_dense : integral_constant<Bool,
	is_dense<T>() && is_real<typename T::value_type>()> {};

template <class T> struct is_comp_dense : integral_constant<Bool,
	is_dense<T>() && is_fp_comp<typename T::value_type>()> {};

// check if is sparse vector/matrix
template <class T> struct is_sparse : integral_constant<Bool,
	is_Diag<T>() || is_MatCoo<T>() || is_MatCooH<T>()> {};

// check if is any slisc container
template <class T> struct is_contain : integral_constant<Bool,
	is_dense<T>() || is_sparse<T>()> {};

template <class T> struct is_real_contain : integral_constant<Bool,
	is_contain<T>() && is_real<typename T::value_type>()> {};

template <class T> struct is_comp_contain : integral_constant<Bool,
	is_contain<T>() && is_fp_comp<typename T::value_type>()> {};

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
template <class T, Int val, SLS_IF(is_arithmetic<typename rm_complex<T>::type>())>
struct Aconst
{
	static constexpr typename rm_complex<T>::type value = val;
};

// Cconst<T, real, imag>::value, where T must be std::complex<>, is a constexpr variable of complex<>
// will not instanciate if T is illegal
template <class T, Int real, Int imag = 0, SLS_IF(is_fp_comp<T>())>
struct Cconst
{
	static constexpr T value = T(real, imag);
};

// === type mapping ===

// implementation of promo_type<T,U>
template <class T1, class T2,
	SLS_IF(is_real<T1>() && is_real<T2>() && type_num<T1>() >= type_num<T2>())>
auto promo_type_fun() { return T1(); }

template <class T1, class T2,
	SLS_IF(is_real<T1>() && is_real<T2>() && type_num<T1>() < type_num<T2>())>
auto promo_type_fun() { return T2(); }

template <class T1, class T2,
	SLS_IF(is_fp_comp<T1>() && is_fp_comp<T2>() && type_num<T1>() >= type_num<T2>())>
auto promo_type_fun() { return T1(); }

template <class T1, class T2,
	SLS_IF(is_fp_comp<T1>() && is_fp_comp<T2>() && type_num<T1>() < type_num<T2>())>
auto promo_type_fun() { return T2(); }

template <class Tr, class Tc,
	SLS_IF(is_real<Tr>() && is_fp_comp<Tc>() && type_num<Tc>() - type_num<Tr>() >= 20)>
auto promo_type_fun() { return Tc(); }

template <class Tc, class Tr,
	SLS_IF(is_real<Tr>() && is_fp_comp<Tc>() && type_num<Tc>() - type_num<Tr>() >= 20)>
auto promo_type_fun() { return Tc(); }

template <class Tr, class Tc,
	SLS_IF(is_real<Tr>() && is_fp_comp<Tc>() && type_num<Tc>() - type_num<Tr>() < 20)>
auto promo_type_fun() { return complex<Tr>(); }

template <class Tc, class Tr,
	SLS_IF(is_real<Tr>() && is_fp_comp<Tc>() && type_num<Tc>() - type_num<Tr>() < 20)>
auto promo_type_fun() { return complex<Tr>(); }

// promo_type<T,U>::type is the smallest type that T, U can losslessly converted to
// (including from integer to floating point conversion)
// will be used as the return type of operator+-*/(T t, U u), etc.
// e.g. Bool + Doub = Doub; Int + Comp = Comp; Doub + complex<float> = Comp;
// see promo_type0() for detail
template <class T1, class T2> struct promo_type
{ using type = decltype(promo_type_fun<T1, T2>()); };

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
