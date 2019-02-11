#pragma once
#include "../SLISC/scalar_arith.h"
#include "../SLISC/matrix.h"
#include "../SLISC/cmat.h"
#include "../SLISC/mat3d.h"
#include "../SLISC/fixsize.h"
#include "../SLISC/sparse.h"

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

	// is_integral
	if (!is_integral<Char>()) error("failed!");
	if (!is_integral<Uchar>()) error("failed!");
	if (!is_integral<Int>()) error("failed!");
	if (!is_integral<Long>()) error("failed!");
	if (!is_integral<Llong>()) error("failed!");
	if (is_integral<float>()) error("failed!");
	if (is_integral<Doub>()) error("failed!");
	if (is_integral<Ldoub>()) error("failed!");

	// is_floating_point
	if (!is_floating_point<float>()) error("failed!");
	if (!is_floating_point<Doub>()) error("failed!");
	if (!is_floating_point<Ldoub>()) error("failed!");
	if (is_floating_point<Char>()) error("failed!");
	if (is_floating_point<Uchar>()) error("failed!");
	if (is_floating_point<Int>()) error("failed!");
	if (is_floating_point<Long>()) error("failed!");
	if (is_floating_point<Llong>()) error("failed!");

	// is_arithmetic
	if (!is_arithmetic<Char>()) error("failed!");
	if (!is_arithmetic<Uchar>()) error("failed!");
	if (!is_arithmetic<Int>()) error("failed!");
	if (!is_arithmetic<Long>()) error("failed!");
	if (!is_arithmetic<float>()) error("failed!");
	if (!is_arithmetic<Doub>()) error("failed!");
	if (!is_arithmetic<Ldoub>()) error("failed!");

	// is_fundamental
	if (!is_fundamental<Char>()) error("failed!");
	if (!is_fundamental<Uchar>()) error("failed!");
	if (!is_fundamental<Int>()) error("failed!");
	if (!is_fundamental<Uint>()) error("failed!");
	if (!is_fundamental<Long>()) error("failed!");
	if (!is_fundamental<Doub>()) error("failed!");
	if (!is_fundamental<Ldoub>()) error("failed!");
	if (is_fundamental<Comp>()) error("failed!");

	// is_signed
	if (!is_signed<Char>()) error("failed!");
	if (is_signed<Uchar>()) error("failed!");
	if (!is_signed<Int>()) error("failed!");
	if (is_signed<Uint>()) error("failed!");
	if (!is_signed<Long>()) error("failed!");
	if (!is_signed<Doub>()) error("failed!");
	if (!is_signed<Ldoub>()) error("failed!");
	if (is_signed<Comp>()) error("failed!");

	// typenum
	if constexpr (type_num<Bool>() != 0) error("failed!");
	if constexpr (type_num<Int>() != 2) error("failed!");
	if constexpr (type_num<Doub>() != 21) error("failed!");
	if constexpr (type_num<Comp>() != 41) error("failed!");
	if constexpr (type_num<std::complex<Uchar>>() != -1) error("failed!");

	// is_comp
	if (!is_comp<Fcomp>()) error("failed!");
	if (!is_comp<Comp>()) error("failed!");
	if (!is_comp<Lcomp>()) error("failed!");
	if (is_comp<complex<Int>>()) error("failed!");
	if (is_comp<Char>()) error("failed!");
	if (is_comp<Int>()) error("failed!");
	if (is_comp<Long>()) error("failed!");
	if (is_comp<Float>()) error("failed!");
	if (is_comp<Doub>()) error("failed!");
	if (is_comp<Ldoub>()) error("failed!");

	// val_type
	if (!is_same<contain_type<VecComp>, Comp>()) error("failed!");
	if (!is_same<contain_type<CmatChar>, Char>()) error("failed!");
	if (!is_same<contain_type<FixVec<Fcomp, 3>>, Fcomp>()) error("failed!");

	// is_contain
	if (!is_contain<Vector<Int>>()) error("failed!");
	if (!is_contain<Matrix<Float>>()) error("failed!");
	if (!is_contain<Cmat<Doub>>()) error("failed!");
	if (!is_contain<FixVec<Comp, 2>>()) error("failed!");
	if (!is_contain<FixCmat<Fcomp, 2, 2>>()) error("failed!");
	if (!is_contain<MatCoo<Char>>()) error("failed!");
	if (!is_contain<MatCooH<Bool>>()) error("failed!");
	if (is_contain<Bool>()) error("failed!");
	if (is_contain<Int>()) error("failed!");
	if (is_contain<Comp>()) error("failed!");
	if (is_contain<vector<Long>>()) error("failed!");

	// is_real_contain
	if (!is_real_contain<Vector<Int>>()) error("failed!");
	if (!is_real_contain<Matrix<Float>>()) error("failed!");
	if (!is_real_contain<Cmat<Doub>>()) error("failed!");
	if (is_real_contain<FixVec<Comp, 2>>()) error("failed!");
	if (is_real_contain<FixCmat<Fcomp, 2, 2>>()) error("failed!");
	if (!is_real_contain<MatCoo<Char>>()) error("failed!");
	if (!is_real_contain<MatCooH<Bool>>()) error("failed!");
	if (is_real_contain<Bool>()) error("failed!");
	if (is_real_contain<Int>()) error("failed!");
	if (is_real_contain<Comp>()) error("failed!");
	if (is_real_contain<vector<Long>>()) error("failed!");

	// is_same_contain
	if constexpr (!is_same_contain<VecComp, VecDoub>()) error("failed!");
	if constexpr (is_same_contain<Comp, VecDoub>()) error("failed!");
	if constexpr (is_same_contain<VecDoub, Float>()) error("failed!");
	if constexpr (is_same_contain<Doub, Float>()) error("failed!");

	// rm_comp
	if (!is_same<rm_comp<Bool>, Bool>()) error("failed!");
	if (!is_same<rm_comp<Int>, Int>()) error("failed!");
	if (!is_same<rm_comp<Char>, Char>()) error("failed!");
	if (!is_same<rm_comp<Float>, Float>()) error("failed!");
	if (!is_same<rm_comp<Doub>, Doub>()) error("failed!");
	if (!is_same<rm_comp<Fcomp>, Float>()) error("failed!");
	if (!is_same<rm_comp<Comp>, Doub>()) error("failed!");
	if (!is_same<rm_comp<Lcomp>, Ldoub>()) error("failed!");

	// Sconst
	if (!is_equiv(Sconst<Int, 0>::value, 0)) error("failed!");
	if (!is_equiv(Sconst<Int, -1>::value, -1)) error("failed!");
	if (!is_equiv(Sconst<Doub, 0>::value, 0.)) error("failed!");
	if (!is_equiv(Sconst<Comp, -1>::value, Comp(-1.))) error("failed!");

	// Rconst
	if (!is_equiv(Rconst<Int, 0>::value, 0)) error("failed!");
	if (!is_equiv(Rconst<Int, -1>::value, -1)) error("failed!");
	if (!is_equiv(Rconst<Doub, 0>::value, 0.)) error("failed!");
	if (!is_equiv(Rconst<Comp, -1>::value, -1.)) error("failed!");

	// Cconst
	if (!is_equiv(Cconst<Fcomp, 0, 1>::value, Fcomp(0, 1))) error("failed!");
	if (!is_equiv(Cconst<Comp, -1>::value, Comp(-1))) error("failed!");
	if (!is_equiv(Cconst<Lcomp, 3, 4>::value, Lcomp(3, 4))) error("failed!");

	// is_promo
	if (!is_promo<Char, Char>()) error("failed!");
	if (!is_promo<Int, Int>()) error("failed!");
	if (!is_promo<Long, Long>()) error("failed!");
	if (!is_promo<Float, Float>()) error("failed!");
	if (!is_promo<Doub, Doub>()) error("failed!");
	if (!is_promo<Comp, Comp>()) error("failed!");
	if (!is_promo<Lcomp, Lcomp>()) error("failed!");

	if (!is_promo<Int, Char>()) error("failed!");
	if (!is_promo<Llong, Int>()) error("failed!");
	if (!is_promo<Float, Llong>()) error("failed!");
	if (!is_promo<Doub, Float>()) error("failed!");
	if (!is_promo<Ldoub, Doub>()) error("failed!");
	if (!is_promo<Comp, Fcomp>()) error("failed!");
	if (!is_promo<Lcomp, Comp>()) error("failed!");

	if (!is_promo<Fcomp, Float>()) error("failed!");
	if (!is_promo<Comp, Float>()) error("failed!");
	if (!is_promo<Lcomp, Float>()) error("failed!");
	if (is_promo<Fcomp, Doub>()) error("failed!");
	if (!is_promo<Comp, Doub>()) error("failed!");
	if (!is_promo<Lcomp, Doub>()) error("failed!");
	if (is_promo<Fcomp, Ldoub>()) error("failed!");
	if (is_promo<Comp, Ldoub>()) error("failed!");
	if (!is_promo<Lcomp, Ldoub>()) error("failed!");

	// promo_type
	if (!is_same<promo_type<Bool, Char>, Char>()) error("failed!");
	if (!is_same<promo_type<Char, Bool>, Char>()) error("failed!");
	if (!is_same<promo_type<Char, Int>, Int>()) error("failed!");
	if (!is_same<promo_type<Int, Char>, Int>()) error("failed!");
	if (!is_same<promo_type<Int, Long>, Long>()) error("failed!");
	if (!is_same<promo_type<Long, Int>, Long>()) error("failed!");
	if (!is_same<promo_type<Long, Float>, Float>()) error("failed!");
	if (!is_same<promo_type<Float, Long>, Float>()) error("failed!");
	if (!is_same<promo_type<Float, Doub>, Doub>()) error("failed!");
	if (!is_same<promo_type<Doub, Float>, Doub>()) error("failed!");
	if (!is_same<promo_type<Doub, Ldoub>, Ldoub>()) error("failed!");
	if (!is_same<promo_type<Ldoub, Doub>, Ldoub>()) error("failed!");
	if (!is_same<promo_type<Doub, Comp>, Comp>()) error("failed!");
	if (!is_same<promo_type<Comp, Doub>, Comp>()) error("failed!");
	if (!is_same<promo_type<Fcomp, Comp>, Comp>()) error("failed!");
	if (!is_same<promo_type<Comp, Fcomp>, Comp>()) error("failed!");
	if (!is_same<promo_type<Doub, Fcomp>, Comp>()) error("failed!");
	if (!is_same<promo_type<Fcomp, Doub>, Comp>()) error("failed!");
	if (!is_same<promo_type<Ldoub, Comp>, Lcomp>()) error("failed!");
	if (!is_same<promo_type<Comp, Ldoub>, Lcomp>()) error("failed!");
	if (!is_same<promo_type<Ldoub, Fcomp>, Lcomp>()) error("failed!");
	if (!is_same<promo_type<Fcomp, Ldoub>, Lcomp>()) error("failed!");

	// num_type
	if (!is_same<num_type<Bool>, Int>()) error("failed!");
	if (!is_same<num_type<Char>, Int>()) error("failed!");
	if (!is_same<num_type<Uchar>, Int>()) error("failed!");
	if (!is_same<num_type<Int>, Int>()) error("failed!");
	if (!is_same<num_type<Long>, Long>()) error("failed!");
	if (!is_same<num_type<Float>, Float>()) error("failed!");
	if (!is_same<num_type<Doub>, Doub>()) error("failed!");
	if (!is_same<num_type<Comp>, Comp>()) error("failed!");
}
