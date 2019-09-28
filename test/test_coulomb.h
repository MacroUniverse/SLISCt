#pragma once
#include "../SLISC/arithmetic.h"
#include "../SLISC/coulomb.h"
#include "../SLISC/time.h"
#include "../SLISC/disp.h"

void test_coulomb()
{
    using std::cout; using std::endl; using std::conj;
    using namespace slisc;
    Int l = 0;
    Doub k = 2.;
    Timer time;

    Long N = 100;
    VecDoub r(N), F(N); linspace(r, 0., 10.);
    //r(0) = 5.671342685370742;
    //r(1) = 5.691382765531062;

    time.tic();
    coulombF(F, l, k, r);
    cout << "time/eval = " << time.toc()/N << endl;
    //disp(r, 15);
    //disp(F, 15);

    // coulomb phase shift
#ifdef SLS_USE_GSL
    Doub ret = coulomb_sigma(3, -2./5);
    if (abs(ret + 0.503297642943251313) > 1e-15)
        SLS_ERR("failed!");
    ret = coulomb_sigma(4, -2. / 9);
    if (abs(ret + 0.3347819876751476) > 1e-15)
        SLS_ERR("failed!");
#endif
}
