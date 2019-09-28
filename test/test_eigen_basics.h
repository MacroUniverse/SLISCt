#pragma once
#include "../SLISC/eigen/eigen_basics.h"

void test_eigen_basics()
{
    using namespace slisc;
    // test mul(MatDoub_O, MatDoub_I, MatDoub_I);
    {
        Doub a_[]{ 1, 2, 3, 4 };
        MatDoub a(2, 2, a_);
        Doub b_[]{ 2, 3, 4, 5 };
        MatDoub b(2, 2, b_);
        Doub c1_[]{ 10, 13, 22, 29 };
        MatDoub c, c1(2, 2, c1_);
        mul(c, a, b);
        if (c != c1) SLS_ERR("failed!");
    }

    // test mul(MatComp_O, MatComp_I, MatComp_I);
    {
        Comp s(1. / sqrt(2.), 1. / sqrt(2.));
        MatComp a(2, 2); a(0) = 1.*s; a(1) = 2.*s; a(2) = 3.*s; a(3) = 4.*s;
        MatComp b(2, 2); b(0) = 2.*s; b(1) = 3.*s; b(2) = 4.*s; b(3) = 5.*s;
        MatComp c, c1(2, 2); c1(0) = 10.*I; c1(1) = 13.*I; c1(2) = 22.*I; c1(3) = 29.*I;
        mul(c, a, b);
        c -= c1;
        if (max_abs(c) > 1e-14)  SLS_ERR("failed!");
    }
}
