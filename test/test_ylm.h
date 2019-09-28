#pragma once
#include "../SLISC/ylm.h"

inline void test_ylm()
{
    using namespace slisc;
    Comp ret = ylm(5, 3, 1., 1);
    if (abs(ret - Comp(0.33207247946604928, -0.047335783997989426)) > 1e-14)
        SLS_ERR("failed");

    ret = ylm(5, -3, 1.7, -1);
    if (abs(ret - Comp(0.284088259485476, -0.040495799317643)) > 1e-14)
        SLS_ERR("failed");

    ret = yl1l2LM(1, 2, 3, 1, 1.1, 2.2, 1.2, 2.3);
    if (abs(ret - Comp(-0.01344167979466, 0.016624624728563)) > 1e-8)
        SLS_ERR("failed");
}
