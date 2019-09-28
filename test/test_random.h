#pragma once
#include "../SLISC/random.h"
#include "../SLISC/arithmetic.h"

void test_random()
{
    using namespace slisc;
    using internal::Ran;
    // same seed
    Ran rangen1(1234), rangen2(1234);
    Long i, N = 10;
    for (i = 0; i < 5; ++i)
        if (rangen1.doub() != rangen2.doub()) SLS_ERR("failed!");
    VecDoub vx(N), vx1(N);
    rand(vx); rand(vx1);
    if (vx == vx1) SLS_ERR("failed!");
    if (max(vx) > 1.) SLS_ERR("failed!");
    if (max(vx) < 0.) SLS_ERR("failed!");
    VecComp vc(N), vc1(N);
    rand(vc); rand(vc1);
    if (vc == vc1) SLS_ERR("failed!");
}
