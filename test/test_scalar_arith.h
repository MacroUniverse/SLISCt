#pragma once
#include "../SLISC/scalar_arith.h"

void test_scalar_arith()
{
    using namespace slisc;

    // is_equiv(s, s)    
    if (!is_equiv(3, 3)) SLS_ERR("failed!");
    if (!is_equiv(3.1, 3.1)) SLS_ERR("failed!");
    { Doub s = 3.1; const Doub &s1 = s;
    if (!is_equiv(s, s1)) SLS_ERR("failed!"); }
    if (is_equiv(3., 3)) SLS_ERR("failed!");
    if (is_equiv(3.1f, 3.1)) SLS_ERR("failed!");

    // to_num(s)
    if (!is_equiv(to_num(false), 0)) SLS_ERR("failed!");
    if (!is_equiv(to_num(true), 1)) SLS_ERR("failed!");
    if (!is_equiv(to_num(Char(-125)), -125)) SLS_ERR("failed!");
    if (!is_equiv(to_num(Uchar(253)), 253)) SLS_ERR("failed!");
    if (!is_equiv(to_num(2.3), 2.3)) SLS_ERR("failed!");
    if (!is_equiv(to_num(Comp(3.1,3.2)), Comp(3.1,3.2))) SLS_ERR("failed!");

    // test scalar functions
    if (!isodd(3) || !isodd(-3) || isodd(4) || isodd(-4))
        SLS_ERR("failed!");
    if (!ispow2(4) || !ispow2(32) || !ispow2(1024) || !ispow2(65536)
        || ispow2(12) || ispow2(48))
        SLS_ERR("failed!");
    if (mod(7, 4) != 3 || mod(6, 4) != 2 || mod(5, 4) != 1 || mod(4, 4) != 0)
        SLS_ERR("failed!");
    if (mod(3, 4) != 3 || mod(2, 4) != 2 || mod(1, 4) != 1 || mod(0, 4) != 0)
        SLS_ERR("failed!");
    if (mod(-1, 4) != 3 || mod(-2, 4) != 2 || mod(-3, 4) != 1 || mod(-4, 4) != 0)
        SLS_ERR("failed!");
    if (mod(-5, 4) != 3 || mod(-6, 4) != 2 || mod(-7, 4) != 1 || mod(-8, 4) != 0)
        SLS_ERR("failed!");
    for (Doub s = -3; s <= 3; s += 0.1) {
        Long n;
        Doub d = 0.5;
        Doub ret = mod(n, s, d);
        if (ret < 0 || abs(ret + n * d - s) > 2e-15)
            SLS_ERR("failed!");
    }

    // operator== involving floating point complex<>
    if (Fcomp(3, 4) != Comp(3, 4)) SLS_ERR("failed!");

    // if (Comp(1,0) == 1.f) SLS_ERR("failed!");

    // operator+,-,*,/ between floatig point complex<> and intrinsic types
    if (!is_equiv(Comp(1.1,2.2) + 3, Comp(4.1, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.1, 2.2) + Char(3), Comp(4.1, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.1, 2.2) + 3.f, Comp(4.1, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.1, 2.2) + 3., Comp(4.1, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(-2 + Comp(1., 2.2), Comp(-1, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Char(-2) + Comp(1., 2.2), Comp(-1, 2.2))) SLS_ERR("failed!");

    if (!is_equiv(Comp(1.1, 2.2) - 3, Comp(-1.9, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.1, 2.2) - 3.f, Comp(-1.9, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.1, 2.2) - 3., Comp(-1.9, 2.2))) SLS_ERR("failed!");
    if (!is_equiv(1 - Comp(1.5, 2.2), Comp(-0.5, -2.2))) SLS_ERR("failed!");
    if (!is_equiv(true - Comp(1.5, 2.2), Comp(-0.5, -2.2))) SLS_ERR("failed!");
    if (!is_equiv(1.f - Comp(1.5, 2.2), Comp(-0.5, -2.2))) SLS_ERR("failed!");
    if (!is_equiv(1. - Comp(1.5, 2.2), Comp(-0.5, -2.2))) SLS_ERR("failed!");

    if (!is_equiv(Comp(1.5, 2.5) * Char(3), Comp(4.5, 7.5))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.5, 2.5) * 3.f, Comp(4.5, 7.5))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.5, 2.5) * 3., Comp(4.5, 7.5))) SLS_ERR("failed!");
    if (!is_equiv(Char(3) * Comp(1.5, 2.5), Comp(4.5, 7.5))) SLS_ERR("failed!");

    if (!is_equiv(Comp(4.5, 7.5) / 3, Comp(1.5, 2.5))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.5, 4.5) / 3.f, Comp(0.5, 1.5))) SLS_ERR("failed!");
    if (!is_equiv(Comp(1.5, 4.5) / 3., Comp(0.5, 1.5))) SLS_ERR("failed!");
    if (!is_equiv(6 / Comp(1, 1), Comp(3, -3))) SLS_ERR("failed!");
    if (!is_equiv(6.f / Comp(1, 1), Comp(3, -3))) SLS_ERR("failed!");
    if (!is_equiv(6. / Comp(1, 1), Comp(3, -3))) SLS_ERR("failed!");
    if (!is_equiv(6. / Fcomp(1, 1), Comp(3, -3))) SLS_ERR("failed!");

    // operator+,-,*,/ between two floatig point std::complex<>    
    if (!is_equiv(Fcomp(1.25, 3.75) + Fcomp(0.75, -1.75), std::complex<float>(2, 2)))
        SLS_ERR("failed!");
    if (!is_equiv(Fcomp(1.25, 3.75) - Comp(-0.75, 0.75), Comp(2, 3)))
        SLS_ERR("failed!");
    if (!is_equiv(Fcomp(1., 2.) / Comp(1., 2.), Comp(1, 0)))
        SLS_ERR("failed!");

    // test factorial()
    if (factorial(5) != 120)
        SLS_ERR("failed!");
    if (factorial(17) != 355687428096000)
        SLS_ERR("failed!");
    if (abs(factorial(100) - 9.33262154439441526e157)/9.33262154439441526e157 > 1e-14)
        SLS_ERR("failed!");
    if (abs(factorial(150) - 5.71338395644585459e262) / 5.71338395644585459e262 > 1e-14)
        SLS_ERR("failed!");
}
