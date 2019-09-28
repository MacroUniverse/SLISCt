#pragma once
#include "../SLISC/disp.h"
#include "../SLISC/mat_fun.h"
#include "../SLISC/slice_arith.h"

void test_mat_fun()
{
    using namespace slisc;
    {
        CmatDoub a(3, 3); // symmetric matrix
        a(0, 0) = a(1, 1) = a(2, 2) = 1.;
        a(0, 1) = a(1, 0) = -2.;
        a(0, 2) = a(2, 0) = 3.;
        a(1, 2) = a(2, 1) = -4.;
        CmatDoub out0(3, 3); // exp(a)
        out0(0, 0) = 3.025725493155477e2;
        out0(1, 1) = 4.036648217771199e2;
        out0(2, 2) = 4.758735878210995e2;
        out0(0, 1) = out0(1, 0) = -3.490811211849774e2;
        out0(0, 2) = out0(2, 0) = 3.792041496895060e2;
        out0(1, 2) = out0(2, 1) = -4.382106846116269e2;

        // test exp_mat_sym();
        {
            CmatDoub out(3, 3);
            exp_mat_sym(out, a, 1);
            out -= out0;
            if (max_abs(out) > 5e-13) SLS_ERR("failed!");
        }

        // test ExpA
        {
            ExpA<Doub> expA(a);
            VecComp v1(3);
            CmatComp out(3, 3); out = 0;
            SvecComp v;
            for (Long j = 0; j < 3; ++j) {
                slice1(v, out, j);
                v1 = 0.; v1[j] = 1;
                expA.mul(v, 1., v1);
            }
            out -= out0;
            if (max_abs(out) > 5e-13) SLS_ERR("failed!");
        }
    }
}
