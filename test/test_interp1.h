#pragma once
#include "../SLISC/interp1.h"

void test_interp1()
{
    using namespace slisc;
    Int i, N0 = 10;
    VecDoub x0(N0); linspace(x0, 0., 2*PI);
    VecDoub y0(x0.size()); sin(y0, x0);
    Spline_interp spline(x0, y0);
    for (i = 0; i < N0; ++i) {
        if (abs(spline.interp(x0[i]) - y0[i]) > 1e-15) SLS_ERR("failed!");;
    }

    linspace(x0, 1., (Doub)N0);
    Spline_interp spline1(x0, x0);
    for (i = 1; i < N0; ++i) {
        if (abs(spline1.interp(x0[i] - 0.5) - (x0[i]-0.5)) > 1e-10) SLS_ERR("failed!");;
    }
}
