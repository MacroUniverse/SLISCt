#pragma once
#include "../SLISC/meta.h"
#include "../SLISC/slisc.h"

inline void test_meta()
{
	using namespace slisc;
	
	// is_same
	if (!is_same<Char, Char>()) error("failed!");
	if (!is_same<Int, Int>()) error("failed!");
	if (!is_same<Long, Long>()) error("failed!");
	if (!is_same<Doub, Doub>()) error("failed!");
	if (is_same<Char, Int>()) error("failed!");
	if (is_same<Int, Doub>()) error("failed!");
	if (is_same<Doub, Comp>()) error("failed!");

	// is_fundamental
	if (!is_fundamental<Char>()) error("failed!");
	if (!is_fundamental<Uchar>()) error("failed!");
	if (!is_fundamental<Int>()) error("failed!");
	if (!is_fundamental<Uint>()) error("failed!");
	if (!is_fundamental<Long>()) error("failed!");
	if (!is_fundamental<Doub>()) error("failed!");
	if (!is_fundamental<Ldoub>()) error("failed!");
	if (is_fundamental<Comp>()) error("failed!");

	// is_floating_point
	if (!is_floating_point<float>()) error("failed!");
	if (!is_floating_point<Doub>()) error("failed!");
	if (!is_floating_point<Ldoub>()) error("failed!");
	if (is_floating_point<Char>()) error("failed!");
	if (is_floating_point<Uchar>()) error("failed!");
	if (is_floating_point<Int>()) error("failed!");
	if (is_floating_point<Long>()) error("failed!");
	if (is_floating_point<Llong>()) error("failed!");

	// is_integral
	if (!is_integral<Char>()) error("failed!");
	if (!is_integral<Uchar>()) error("failed!");
	if (!is_integral<Int>()) error("failed!");
	if (!is_integral<Long>()) error("failed!");
	if (!is_integral<Llong>()) error("failed!");
	if (is_integral<float>()) error("failed!");
	if (is_integral<Doub>()) error("failed!");
	if (is_integral<Ldoub>()) error("failed!");

	// is_complex
	if (!is_complex<std::complex<Int>>()) error("failed!");
	if (!is_complex<std::complex<float>>()) error("failed!");
	if (!is_complex<Comp>()) error("failed!");
	if (!is_complex<std::complex<Ldoub>>()) error("failed!");
	if (is_complex<Char>()) error("failed!");
	if (is_complex<Int>()) error("failed!");
	if (is_complex<Long>()) error("failed!");
	if (is_complex<float>()) error("failed!");
	if (is_complex<Doub>()) error("failed!");
	if (is_complex<Ldoub>()) error("failed!");

	// typenum
	if constexpr (type_num<Bool>() != 1) error("failed!");
	if constexpr (type_num<Int>() != 4) error("failed!");
	if constexpr (type_num<Doub>() != 12) error("failed!");
	if constexpr (type_num<std::complex<Uchar>>() != -1) error("failed!");

	// promo_type
	if (!is_same<promo_type<Bool, Char>::type, Char>()) error("failed!");
	if (!is_same<promo_type<Char, Int>::type, Int>()) error("failed!");
	if (!is_same<promo_type<Int, Long>::type, Long>()) error("failed!");
	if (!is_same<promo_type<Long, float>::type, float>()) error("failed!");
	if (!is_same<promo_type<float, Doub>::type, Doub>()) error("failed!");
	if (!is_same<promo_type<Doub, Ldoub>::type, Ldoub>()) error("failed!");
	if (!is_same<promo_type<Doub, Comp>::type, Comp>()) error("failed!");
	if (!is_same<promo_type<std::complex<float>, Comp>::type, Comp>()) error("failed!");

	// to_num_t
	if (!is_same<to_num_t<Bool>::type, Int>()) error("failed!");
	if (!is_same<to_num_t<Char>::type, Int>()) error("failed!");
	if (!is_same<to_num_t<Uchar>::type, Int>()) error("failed!");
	if (!is_same<to_num_t<Int>::type, Int>()) error("failed!");
	if (!is_same<to_num_t<Long>::type, Long>()) error("failed!");
	if (!is_same<to_num_t<float>::type, float>()) error("failed!");
	if (!is_same<to_num_t<Doub>::type, Doub>()) error("failed!");
	if (!is_same<to_num_t<Comp>::type, Comp>()) error("failed!");

	// rm_comp
	if (!is_same<rm_comp<Bool>::type, Bool>()) error("failed!");
	if (!is_same<rm_comp<Int>::type, Int>()) error("failed!");
	if (!is_same<rm_comp<Char>::type, Char>()) error("failed!");
	if (!is_same<rm_comp<float>::type, float>()) error("failed!");
	if (!is_same<rm_comp<Doub>::type, Doub>()) error("failed!");
	if (!is_same<rm_comp<std::complex<float>>::type, float>()) error("failed!");
	if (!is_same<rm_comp<Comp>::type, Doub>()) error("failed!");
	if (!is_same<rm_comp<std::complex<Ldoub>>::type, Ldoub>()) error("failed!");
}
