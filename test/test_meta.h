#pragma once
#include "../SLISC/scalar_arith.h"
#include "../SLISC/matrix.h"
#include "../SLISC/cmat.h"
#include "../SLISC/mat3d.h"
#include "../SLISC/fixsize.h"
#include "../SLISC/matcooh.h"

// #include "../SLISC/slisc.h"

inline void test_meta()
{
    using namespace slisc;
    
    // is_same
    if (!is_same<Char, Char>()) SLS_ERR("failed!");
    if (!is_same<Int, Int>()) SLS_ERR("failed!");
    if (!is_same<Long, Long>()) SLS_ERR("failed!");
    if (!is_same<Doub, Doub>()) SLS_ERR("failed!");
    if (is_same<Char, Int>()) SLS_ERR("failed!");
    if (is_same<Int, Doub>()) SLS_ERR("failed!");
    if (is_same<Doub, Comp>()) SLS_ERR("failed!");

    // is_integral
    if (!is_integral<Char>()) SLS_ERR("failed!");
    if (!is_integral<Uchar>()) SLS_ERR("failed!");
    if (!is_integral<Int>()) SLS_ERR("failed!");
    if (!is_integral<Long>()) SLS_ERR("failed!");
    if (!is_integral<Llong>()) SLS_ERR("failed!");
    if (is_integral<Float>()) SLS_ERR("failed!");
    if (is_integral<Doub>()) SLS_ERR("failed!");
    if (is_integral<Ldoub>()) SLS_ERR("failed!");

    // is_fpt
    if (!is_fpt<Float>()) SLS_ERR("failed!");
    if (!is_fpt<Doub>()) SLS_ERR("failed!");
    if (!is_fpt<Ldoub>()) SLS_ERR("failed!");
    if (is_fpt<Char>()) SLS_ERR("failed!");
    if (is_fpt<Uchar>()) SLS_ERR("failed!");
    if (is_fpt<Int>()) SLS_ERR("failed!");
    if (is_fpt<Long>()) SLS_ERR("failed!");
    if (is_fpt<Llong>()) SLS_ERR("failed!");

    // is_arithmetic
    if (!is_arithmetic<Char>()) SLS_ERR("failed!");
    if (!is_arithmetic<Uchar>()) SLS_ERR("failed!");
    if (!is_arithmetic<Int>()) SLS_ERR("failed!");
    if (!is_arithmetic<Long>()) SLS_ERR("failed!");
    if (!is_arithmetic<Float>()) SLS_ERR("failed!");
    if (!is_arithmetic<Doub>()) SLS_ERR("failed!");
    if (!is_arithmetic<Ldoub>()) SLS_ERR("failed!");

    // is_fundamental
    if (!is_fundamental<Char>()) SLS_ERR("failed!");
    if (!is_fundamental<Uchar>()) SLS_ERR("failed!");
    if (!is_fundamental<Int>()) SLS_ERR("failed!");
    if (!is_fundamental<Uint>()) SLS_ERR("failed!");
    if (!is_fundamental<Long>()) SLS_ERR("failed!");
    if (!is_fundamental<Doub>()) SLS_ERR("failed!");
    if (!is_fundamental<Ldoub>()) SLS_ERR("failed!");
    if (is_fundamental<Comp>()) SLS_ERR("failed!");

    // is_signed
    if (!is_signed<Char>()) SLS_ERR("failed!");
    if (is_signed<Uchar>()) SLS_ERR("failed!");
    if (!is_signed<Int>()) SLS_ERR("failed!");
    if (is_signed<Uint>()) SLS_ERR("failed!");
    if (!is_signed<Long>()) SLS_ERR("failed!");
    if (!is_signed<Doub>()) SLS_ERR("failed!");
    if (!is_signed<Ldoub>()) SLS_ERR("failed!");
    if (is_signed<Comp>()) SLS_ERR("failed!");

    // typenum
    if (type_num<Bool>() != 0) SLS_ERR("failed!");
    if (type_num<Int>() != 2) SLS_ERR("failed!");
    if (type_num<Doub>() != 21) SLS_ERR("failed!");
    if (type_num<Comp>() != 41) SLS_ERR("failed!");
    if (type_num<std::complex<Uchar>>() != -1) SLS_ERR("failed!");

    // is_comp
    if (!is_comp<Fcomp>()) SLS_ERR("failed!");
    if (!is_comp<Comp>()) SLS_ERR("failed!");
    if (!is_comp<Lcomp>()) SLS_ERR("failed!");
    if (is_comp<complex<Int>>()) SLS_ERR("failed!");
    if (is_comp<Char>()) SLS_ERR("failed!");
    if (is_comp<Int>()) SLS_ERR("failed!");
    if (is_comp<Long>()) SLS_ERR("failed!");
    if (is_comp<Float>()) SLS_ERR("failed!");
    if (is_comp<Doub>()) SLS_ERR("failed!");
    if (is_comp<Ldoub>()) SLS_ERR("failed!");

    // val_type
    if (!is_same<contain_type<VecComp>, Comp>()) SLS_ERR("failed!");
    if (!is_same<contain_type<CmatChar>, Char>()) SLS_ERR("failed!");
    if (!is_same<contain_type<FixVec<Fcomp, 3>>, Fcomp>()) SLS_ERR("failed!");

    // is_contain
    if (!is_contain<Vector<Int>>()) SLS_ERR("failed!");
    if (!is_contain<Matrix<Float>>()) SLS_ERR("failed!");
    if (!is_contain<Cmat<Doub>>()) SLS_ERR("failed!");
    if (!is_contain<FixVec<Comp, 2>>()) SLS_ERR("failed!");
    if (!is_contain<FixCmat<Fcomp, 2, 2>>()) SLS_ERR("failed!");
    if (!is_contain<MatCoo<Char>>()) SLS_ERR("failed!");
    if (!is_contain<MatCooH<Bool>>()) SLS_ERR("failed!");
    if (is_contain<Bool>()) SLS_ERR("failed!");
    if (is_contain<Int>()) SLS_ERR("failed!");
    if (is_contain<Comp>()) SLS_ERR("failed!");
    if (is_contain<vector<Long>>()) SLS_ERR("failed!");

    // is_real_contain
    if (!is_real_contain<Vector<Int>>()) SLS_ERR("failed!");
    if (!is_real_contain<Matrix<Float>>()) SLS_ERR("failed!");
    if (!is_real_contain<Cmat<Doub>>()) SLS_ERR("failed!");
    if (is_real_contain<FixVec<Comp, 2>>()) SLS_ERR("failed!");
    if (is_real_contain<FixCmat<Fcomp, 2, 2>>()) SLS_ERR("failed!");
    if (!is_real_contain<MatCoo<Char>>()) SLS_ERR("failed!");
    if (!is_real_contain<MatCooH<Bool>>()) SLS_ERR("failed!");
    if (is_real_contain<Bool>()) SLS_ERR("failed!");
    if (is_real_contain<Int>()) SLS_ERR("failed!");
    if (is_real_contain<Comp>()) SLS_ERR("failed!");
    if (is_real_contain<vector<Long>>()) SLS_ERR("failed!");

    // is_same_contain
    if (!is_same_contain<VecComp, VecDoub>()) SLS_ERR("failed!");
    if (is_same_contain<Comp, VecDoub>()) SLS_ERR("failed!");
    if (is_same_contain<VecDoub, Float>()) SLS_ERR("failed!");
    if (is_same_contain<Doub, Float>()) SLS_ERR("failed!");

    // rm_comp
    if (!is_same<rm_comp<Bool>, Bool>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Int>, Int>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Char>, Char>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Float>, Float>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Doub>, Doub>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Fcomp>, Float>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Comp>, Doub>()) SLS_ERR("failed!");
    if (!is_same<rm_comp<Lcomp>, Ldoub>()) SLS_ERR("failed!");

    // Sconst
    if (!is_equiv(Sconst<Int, 0>::value, 0)) SLS_ERR("failed!");
    if (!is_equiv(Sconst<Int, -1>::value, -1)) SLS_ERR("failed!");
    if (!is_equiv(Sconst<Doub, 0>::value, 0.)) SLS_ERR("failed!");
    if (!is_equiv(Sconst<Comp, -1>::value, Comp(-1.))) SLS_ERR("failed!");

    // Rconst
    if (!is_equiv(Rconst<Int, 0>::value, 0)) SLS_ERR("failed!");
    if (!is_equiv(Rconst<Int, -1>::value, -1)) SLS_ERR("failed!");
    if (!is_equiv(Rconst<Doub, 0>::value, 0.)) SLS_ERR("failed!");
    if (!is_equiv(Rconst<Comp, -1>::value, -1.)) SLS_ERR("failed!");

    // Cconst
    if (!is_equiv(Cconst<Fcomp, 0, 1>::value, Fcomp(0, 1))) SLS_ERR("failed!");
    if (!is_equiv(Cconst<Comp, -1>::value, Comp(-1))) SLS_ERR("failed!");
    if (!is_equiv(Cconst<Lcomp, 3, 4>::value, Lcomp(3, 4))) SLS_ERR("failed!");

    // is_promo
    if (!is_promo<Char, Char>()) SLS_ERR("failed!");
    if (!is_promo<Int, Int>()) SLS_ERR("failed!");
    if (!is_promo<Long, Long>()) SLS_ERR("failed!");
    if (!is_promo<Float, Float>()) SLS_ERR("failed!");
    if (!is_promo<Doub, Doub>()) SLS_ERR("failed!");
    if (!is_promo<Comp, Comp>()) SLS_ERR("failed!");
    if (!is_promo<Lcomp, Lcomp>()) SLS_ERR("failed!");

    if (!is_promo<Int, Char>()) SLS_ERR("failed!");
    if (!is_promo<Llong, Int>()) SLS_ERR("failed!");
    if (!is_promo<Float, Llong>()) SLS_ERR("failed!");
    if (!is_promo<Doub, Float>()) SLS_ERR("failed!");
    if (!is_promo<Ldoub, Doub>()) SLS_ERR("failed!");
    if (!is_promo<Comp, Fcomp>()) SLS_ERR("failed!");
    if (!is_promo<Lcomp, Comp>()) SLS_ERR("failed!");

    if (!is_promo<Fcomp, Float>()) SLS_ERR("failed!");
    if (!is_promo<Comp, Float>()) SLS_ERR("failed!");
    if (!is_promo<Lcomp, Float>()) SLS_ERR("failed!");
    if (is_promo<Fcomp, Doub>()) SLS_ERR("failed!");
    if (!is_promo<Comp, Doub>()) SLS_ERR("failed!");
    if (!is_promo<Lcomp, Doub>()) SLS_ERR("failed!");
    if (is_promo<Fcomp, Ldoub>()) SLS_ERR("failed!");
    if (is_promo<Comp, Ldoub>()) SLS_ERR("failed!");
    if (!is_promo<Lcomp, Ldoub>()) SLS_ERR("failed!");

    // promo_type
    if (!is_same<promo_type<Bool, Char>, Char>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Char, Bool>, Char>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Char, Int>, Int>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Int, Char>, Int>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Int, Long>, Long>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Long, Int>, Long>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Long, Float>, Float>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Float, Long>, Float>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Float, Doub>, Doub>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Doub, Float>, Doub>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Doub, Ldoub>, Ldoub>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Ldoub, Doub>, Ldoub>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Doub, Comp>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Comp, Doub>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Fcomp, Comp>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Comp, Fcomp>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Doub, Fcomp>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Fcomp, Doub>, Comp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Ldoub, Comp>, Lcomp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Comp, Ldoub>, Lcomp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Ldoub, Fcomp>, Lcomp>()) SLS_ERR("failed!");
    if (!is_same<promo_type<Fcomp, Ldoub>, Lcomp>()) SLS_ERR("failed!");

    // num_type
    if (!is_same<num_type<Bool>, Int>()) SLS_ERR("failed!");
    if (!is_same<num_type<Char>, Int>()) SLS_ERR("failed!");
    if (!is_same<num_type<Uchar>, Int>()) SLS_ERR("failed!");
    if (!is_same<num_type<Int>, Int>()) SLS_ERR("failed!");
    if (!is_same<num_type<Long>, Long>()) SLS_ERR("failed!");
    if (!is_same<num_type<Float>, Float>()) SLS_ERR("failed!");
    if (!is_same<num_type<Doub>, Doub>()) SLS_ERR("failed!");
    if (!is_same<num_type<Comp>, Comp>()) SLS_ERR("failed!");
}
