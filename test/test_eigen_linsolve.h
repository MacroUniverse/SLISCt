#pragma once
#include "../SLISC/eigen/eigen_linsolve.h"

// using eigen_linsolve.h will cause gcc error????
void test_eigen_linsolve()
{
    using namespace slisc;
#ifdef _MSC_VER
    {
        MatDoub A(2, 2), B(2, 2);
        A(0) = 1.; A(1) = 2.; A(2) = 3.; A(3) = 4.;
        plus(B, A, 1.);
        VecDoub x(2), x1(2), x2(2);
        x1[0] = 1.; x1[1] = 2.;
        plus(x2, x1, 1.);
        VecDoub y1(2), y2(2);
        mul(y1, A, x1); mul(y2, A, x2);
        HouseholderQrDoub qr(A);
        qr.solve(x, y1);
        x -= x1; abs(x);
        if (max(x) > 2e-14) SLS_ERR("failed!");
        qr.solve(x, y2);
        x -= x2; abs(x);
        if (max(x) > 2e-14) SLS_ERR("failed!");
        mul(y1, B, x1);
        qr.compute(B);
        qr.solve(x, y1);
        x -= x1; abs(x);
        if (max(x) > 2.5e-14) SLS_ERR("failed!");
    }

    {
        MatComp A(2, 2), B(2, 2);
        A(0) = Comp(1., 2.); A(1) = Comp(2., 3.); A(2) = Comp(3., 4.); A(3) = Comp(4., 5.);
        plus(B, A, 1.);
        VecComp x(2), x1(2), x2(2);
        x1[0] = Comp(1., 2.); x1[1] = Comp(3., 4.);
        plus(x2, x1, 1.);
        VecComp y1(2), y2(2);
        mul(y1, A, x1); mul(y2, A, x2);
        HouseholderQrComp qr(A);
        qr.solve(x, y1);
        x -= x1;
        if (max_abs(x) > 2e-14) SLS_ERR("failed!");
        qr.solve(x, y2);
        x -= x2;
        if (max_abs(x) > 2e-14) SLS_ERR("failed!");
        mul(y1, B, x1);
        qr.compute(B);
        qr.solve(x, y1);
        x -= x1;
        if (max_abs(x) > 2.5e-14) SLS_ERR("failed!");
    }

#else
    std::cout << "Not in Visual Studio, not testing eigen_linsolve!" << std::endl;
#endif
}
