#pragma once
#include "../SLISC/global.h"
#include "../SLISC/meta.h"
// #include "../SLISC/slisc.h"

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

	// is_fp_complex
	if (!is_fp_comp<std::complex<float>>()) error("failed!");
	if (!is_fp_comp<Comp>()) error("failed!");
	if (!is_fp_comp<std::complex<Ldoub>>()) error("failed!");
	if (is_fp_comp<std::complex<Int>>()) error("failed!");
	if (is_fp_comp<Char>()) error("failed!");
	if (is_fp_comp<Int>()) error("failed!");
	if (is_fp_comp<Long>()) error("failed!");
	if (is_fp_comp<float>()) error("failed!");
	if (is_fp_comp<Doub>()) error("failed!");
	if (is_fp_comp<Ldoub>()) error("failed!");

	// is_contain
	if (!is_contain<Vector<Int>>()) error("failed!");
	if (!is_contain<Matrix<float>>()) error("failed!");
	if (!is_contain<Cmat<Doub>>()) error("failed!");
	if (!is_contain<FixVec<Comp, 2>>()) error("failed!");
	if (!is_contain<FixCmat<std::complex<float>, 2, 2>>()) error("failed!");
	if (!is_contain<MatCoo<Char>>()) error("failed!");
	if (!is_contain<MatCooH<Bool>>()) error("failed!");
	if (is_contain<Bool>()) error("failed!");
	if (is_contain<Int>()) error("failed!");
	if (is_contain<Comp>()) error("failed!");
	if (is_contain<std::vector<Long>>()) error("failed!");

	// typenum
	if constexpr (type_num<Bool>() != 0) error("failed!");
	if constexpr (type_num<Int>() != 2) error("failed!");
	if constexpr (type_num<Doub>() != 21) error("failed!");
	if constexpr (type_num<Comp>() != 41) error("failed!");
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
	if (!is_same<num_type<Bool>::type, Int>()) error("failed!");
	if (!is_same<num_type<Char>::type, Int>()) error("failed!");
	if (!is_same<num_type<Uchar>::type, Int>()) error("failed!");
	if (!is_same<num_type<Int>::type, Int>()) error("failed!");
	if (!is_same<num_type<Long>::type, Long>()) error("failed!");
	if (!is_same<num_type<float>::type, float>()) error("failed!");
	if (!is_same<num_type<Doub>::type, Doub>()) error("failed!");
	if (!is_same<num_type<Comp>::type, Comp>()) error("failed!");

	// rm_comp
	if (!is_same<rm_complex<Bool>::type, Bool>()) error("failed!");
	if (!is_same<rm_complex<Int>::type, Int>()) error("failed!");
	if (!is_same<rm_complex<Char>::type, Char>()) error("failed!");
	if (!is_same<rm_complex<float>::type, float>()) error("failed!");
	if (!is_same<rm_complex<Doub>::type, Doub>()) error("failed!");
	if (!is_same<rm_complex<std::complex<float>>::type, float>()) error("failed!");
	if (!is_same<rm_complex<Comp>::type, Doub>()) error("failed!");
	if (!is_same<rm_complex<std::complex<Ldoub>>::type, Ldoub>()) error("failed!");

	// is_same_contain
	if (!is_same_contain<VecComp, VecDoub>()) error("failed!");
}
