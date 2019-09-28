#pragma once
#include "../SLISC/eigen/eigen_fft.h"

void test_eigen_fft()
{
    using namespace slisc;
    VecComp f0(4); f0[0] = 1.; f0[1] = I; f0[2] = -1.; f0[3] = -I;
    const VecComp f(4, f0.ptr());
    VecComp g(4), f1(4);
    FFT fft;
    fft.fwd(g, f);
    f1[0] = 0.; f1[1] = 4.; f1[2] = 0.; f1[3] = 0.;
    g -= f1;
    if (max_abs(g) > 1e-15) SLS_ERR("failed!");
    fft.fwd(g, f);
    fft.inv(f1, g);
    f1 /= 4.; f1 -= f0;
    if (max_abs(f1) > 1e-15) SLS_ERR("failed!");
}
