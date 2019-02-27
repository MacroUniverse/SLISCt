// meta-programming utilities

#pragma once
#include "global.h"
#include <type_traits>

// using variable argument macro to allow parsing of ",".
// otherwise, "," will separate a single argument into multiple arguments
#define SLS_IF_HELPER(cond) typename std::enable_if<(bool)(cond), Int>::type = 0
#define SLS_IF(...) SLS_IF_HELPER(__VA_ARGS__)

namespace slisc {

// type traits

// true_type::value is true, false_type::value is false
// true_type() and false_type() can be implicitly converted to true and false
using std::true_type; using std::false_type;

// integral_constant<Bool, true> is true_type
// integral_constant<Bool, false> is false_type
using std::integral_constant;

// is_same<U,T> checks if U = T
template <class T1, class T2>
constexpr Bool is_same()
{
	return std::is_same<T1, T2>::value;
}

// is_integral<T> checks if T = bool, character types, integer types
template <class T>
constexpr Bool is_integral()
{
	return std::is_integral<T>::value;
}

// is_floating_point<T> checks if T = float, double, long double
template <class T>
constexpr Bool is_floating_point()
{
	return std::is_floating_point<T>::value;
}

// is_arithmetic<T> checks if T = integral type || floating point type
template <class T>
constexpr Bool is_arithmetic()
{
	return std::is_arithmetic<T>::value;
}

// is_fundamental<T> checks if T = arithmetic type || void || null pointer
template <class T>
constexpr Bool is_fundamental()
{
	return std::is_fundamental<T>::value;
}

// is_signed<T> checks if T(-1) < T(0)
template <class T>
constexpr Bool is_signed()
{
	return std::is_signed<T>::value;
}

// check if T is a specific type

template<class T>
constexpr Bool is_Bool()
{
	return is_same<T, Bool>();
}

template<class T>
constexpr Bool is_Char()
{
	return is_same<T, Char>();
}

template<class T>
constexpr Bool is_Int()
{
	return is_same<T, Int>();
}

template<class T>
constexpr Bool is_Llong()
{
	return is_same<T, Llong>();
}

template<class T>
constexpr Bool is_Float()
{
	return is_same<T, Float>();
}

template<class T>
constexpr Bool is_Doub()
{
	return is_same<T, Doub>();
}

template<class T>
constexpr Bool is_Ldoub()
{
	return is_same<T, Ldoub>();
}

template<class T>
constexpr Bool is_Fcomp()
{
	return is_same<T, Fcomp>();
}

template<class T>
constexpr Bool is_Comp()
{
	return is_same<T, Comp>();
}

template<class T>
constexpr Bool is_Lcomp()
{
	return is_same<T, Lcomp>();
}

// type_num<T>() maps each scalr type to a unique number
// a scalar type is defined as having a non-negative type number
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
	else {
		if (is_Fcomp<T>()) return 40;
		if (is_Comp<T>()) return 41;
		if (is_Lcomp<T>()) return 42;
	}
	return -1;
}

// is_scalar<T> checks if is a scalar
template<class T>
constexpr Bool is_scalar()
{
	return type_num<T>() >= 0;
}


// is_real<T> checks if is a real scalar
template<class T>
constexpr Bool is_real()
{
	return type_num<T>() >= 0 && type_num<T>() < 40;
}

// is_comp<T> checks if T is a scalar of complex<> type
template<class T>
constexpr Bool is_comp()
{
	return type_num<T>() >= 40;
}

// check if is a specific container type

template <class T> struct is_Vector_imp : false_type {};
template <class T> struct is_Vector_imp<Vector<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_Vector()
{
	return is_Vector_imp<T>();
}

template <class T> struct is_Matrix_imp : false_type {};
template <class T> struct is_Matrix_imp<Matrix<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_Matrix()
{
	return is_Matrix_imp<T>();
}

template <class T> struct is_Cmat_imp : false_type {};
template <class T> struct is_Cmat_imp<Cmat<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_Cmat()
{
	return is_Cmat_imp<T>();
}

template <class T> struct is_FixVec_imp : false_type {};
template <class T, Long N> struct is_FixVec_imp<FixVec<T, N>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_FixVec()
{
	return is_FixVec_imp<T>();
}

template <class T> struct is_FixCmat_imp : false_type {};
template <class T, Long Nr, Long Nc> struct is_FixCmat_imp<FixCmat<T, Nr, Nc>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_FixCmat()
{
	return is_FixCmat_imp<T>();
}

template <class T> struct is_Mat3d_imp : false_type {};
template <class T> struct is_Mat3d_imp<Mat3d<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_Mat3d()
{
	return is_Mat3d_imp<T>();
}

template <class T> struct is_Diag_imp : false_type {};
template <class T> struct is_Diag_imp<Diag<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_Diag()
{
	return is_Diag_imp<T>();
}

template <class T> struct is_MatCoo_imp : false_type {};
template <class T> struct is_MatCoo_imp<MatCoo<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_MatCoo()
{
	return is_MatCoo_imp<T>();
}

template <class T> struct is_MatCooH_imp : false_type {};
template <class T> struct is_MatCooH_imp<MatCooH<T>> : integral_constant<Bool, is_scalar<T>()> {};
template<class T>
constexpr Bool is_MatCooH()
{
	return is_MatCooH_imp<T>();
}

// check if is fixed-size container
template <class T> constexpr Bool is_fixed()
{
	return is_FixVec<T>() || is_FixCmat<T>();
}

// check if is dense matrix (2D)
template <class T> constexpr Bool is_dense_mat()
{
	return is_Matrix<T>() || is_Cmat<T>() || is_FixCmat<T>();
}

// check if is dense container (including fixed-size)
template <class T> constexpr Bool is_dense()
{
	return is_Vector<T>() || is_FixVec<T>() || is_dense_mat<T>() || is_Mat3d<T>();
}

// check if is sparse vector/matrix
template <class T> constexpr Bool is_sparse()
{
	return is_Diag<T>() || is_MatCoo<T>() || is_MatCooH<T>();
}

template <class T>
constexpr Int contain_num()
{
	if (is_Vector<T>()) return 0;
	else if (is_Matrix<T>()) return 1;
	else if (is_Cmat<T>()) return 2;
	else if (is_Mat3d<T>()) return 3;

	else if (is_FixVec<T>()) return 20;
	else if (is_FixCmat<T>()) return 22;

	else if (is_Diag<T>()) return 31;
	else if (is_MatCoo<T>()) return 32;
	else if (is_MatCooH<T>()) return 33;

	return -1;
}

// check if is a slisc container
template <class T> constexpr Bool is_contain()
{
	return contain_num<T>() >= 0;
}

// contain_type<T> is T::value_type, if T is a container
// otherwise, contain_type<T> is T
template <class T, SLS_IF(!is_contain<T>())>
constexpr T contain_type_fun() { return T(); };

template <class T, SLS_IF(is_contain<T>())>
constexpr auto contain_type_fun() { return T::value_type(); };

template <class T> using contain_type = decltype(contain_type_fun<T>());

template <class T> constexpr Bool is_real_dense()
{
	return is_dense<T>() && is_real<contain_type<T>>();
}

template <class T> constexpr Bool is_comp_dense()
{
	return is_dense<T>() && is_comp<contain_type<T>>();
}

// check if is a real slisc container
template <class T> constexpr Bool is_real_contain()
{
	return is_contain<T>() && is_real<contain_type<T>>();
}

// check if is a complex slisc container
template <class T> constexpr Bool is_comp_contain()
{
	return is_contain<T>() && is_comp<contain_type<T>>();
}

// check if two containers are the same (value_type can be different)
template <class T1, class T2> constexpr Bool is_same_contain()
{
	return is_contain<T1>() && contain_num<T1>() == contain_num<T2>();
}

// for Tc = complex<Tr>
// rm_comp<Tc> is Tr
template <class T> struct rm_comp_imp { typedef T type; };
template <class T> struct rm_comp_imp<complex<T>> { typedef T type; };
template <class T> using rm_comp = typename rm_comp_imp<T>::type;

// === get static constexpr values ===

// Sconst<T, val>::value is a constexpr variable of a scalar
template <class T, Int val, SLS_IF(is_scalar<T>())>
struct Sconst
{
	static constexpr T value = val;
};

// Aconst<T, val>::value is a constexpr variable of a real scalar
// T can be any scalar
template <class T, Int val, SLS_IF(is_scalar<T>())>
struct Rconst
{
	static constexpr rm_comp<T> value = val;
};

// Cconst<T, real, imag = 0>::value is a constexpr variable of a complex scalar
template <class T, Int real, Int imag = 0, SLS_IF(is_comp<T>())>
struct Cconst
{
	static constexpr T value = T(real, imag);
};

// === type mapping ===

// is_promo<T1,T2> checks if T2 can be lesslessly converted to T1
// (including T1 = T2)
// (including from integer to floating point conversions)
// might be used to enable operator=,+=,-=,*=,/= etc.
template <class T1, class T2>
constexpr Bool is_promo_fun()
{
	if (is_real<T2>()) {
		if (is_real<T1>()) {
			if (type_num<T1>() >= type_num<T2>())
				return true;
		}
		else { // is_comp<T1>()
			if (type_num<T1>() - type_num<T2>() >= 20)
				return true;
		}
	}
	else { // is_comp<T2>()
		if (is_comp<T1>()) {
			if (type_num<T1>() >= type_num<T2>())
				return true;
		}
	}
	return false;
}

template <class T1, class T2>
constexpr Bool is_promo()
{
	return is_promo_fun<T1, T2>();
}

// promo_type<T1,T2> is the smallest type that both T1, T2 can losslessly converted to
// (including from integer to floating point conversions)
// might be used as the return type of operator+-*/(T1 t1, T2 t2), etc.
// e.g. Bool + Doub = Doub; Int + Comp = Comp; Doub + Fcomp = Comp;

template <class T1, class T2,
	SLS_IF(is_real<T1>() && is_real<T2>() && type_num<T1>() >= type_num<T2>())>
auto promo_type_fun() { return T1(); }

template <class T1, class T2,
	SLS_IF(is_real<T1>() && is_real<T2>() && type_num<T1>() < type_num<T2>())>
auto promo_type_fun() { return T2(); }

template <class T1, class T2,
	SLS_IF(is_comp<T1>() && is_comp<T2>() && type_num<T1>() >= type_num<T2>())>
auto promo_type_fun() { return T1(); }

template <class T1, class T2,
	SLS_IF(is_comp<T1>() && is_comp<T2>() && type_num<T1>() < type_num<T2>())>
auto promo_type_fun() { return T2(); }

template <class Tr, class Tc,
	SLS_IF(is_real<Tr>() && is_comp<Tc>() && type_num<Tc>() - type_num<Tr>() >= 20)>
auto promo_type_fun() { return Tc(); }

template <class Tc, class Tr,
	SLS_IF(is_real<Tr>() && is_comp<Tc>() && type_num<Tc>() - type_num<Tr>() >= 20)>
auto promo_type_fun() { return Tc(); }

template <class Tr, class Tc,
	SLS_IF(is_real<Tr>() && is_comp<Tc>() && type_num<Tc>() - type_num<Tr>() < 20)>
auto promo_type_fun() { return complex<Tr>(); }

template <class Tc, class Tr,
	SLS_IF(is_real<Tr>() && is_comp<Tc>() && type_num<Tc>() - type_num<Tr>() < 20)>
auto promo_type_fun() { return complex<Tr>(); }

template <class T1, class T2> struct promo_type_imp
{ using type = decltype(promo_type_fun<T1, T2>()); };

template <class T1, class T2> using promo_type = typename promo_type_imp<T1, T2>::type;

// to_num_t<T> is Int if T = bool, character types
// otherwise is T
template <typename T> struct num_type_imp { typedef T type; };
template<> struct num_type_imp<Bool> { typedef Int type; };
template<> struct num_type_imp<Char> { typedef Int type; };
template<> struct num_type_imp<Uchar> { typedef Int type; };

template <typename T> using num_type = typename num_type_imp<T>::type;

} // namespace slisc
