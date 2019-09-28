#pragma once
#include "../SLISC/cmatobd.h"

void test_cmatobd()
{
    using namespace slisc;

    Cmat3Int a0(3, 3, 4);
    for (Long k = 0; k < 4; ++k) {
        ScmatInt sli; slice12(sli, a0, k);
        linspace(sli, 1, 9); sli.end() = 1;
    }
    CmobdInt a(3, 4); a = a0;

    // construct from Cmat3d
    {
        VecInt x(7), y(7); linspace(x, 1, 7);
        mul(y, a, x);
        if (y[0] != 21 || y[1] != 48 || y[2] != 51 || y[3] != 90 ||
            y[4] != 81 || y[5] != 76 || y[6] != 47)
            SLS_ERR("failed!");
    }
    
    // construct from Mcoo
    {
        McooInt a1(7, 7, 49);
        for (Long i = 0; i < 7; ++i) {
            for (Long j = 0; j < 7; ++j) {
                if (a(i, j) != 0)
                    a1.push(a(i, j), i, j);
            }
        }
        CmatObd<Int> a2(3, 4); a2 = a1;

        VecInt x(7), y(7); linspace(x, 1, 7);
        mul(y, a2, x);
        if (y[0] != 21 || y[1] != 48 || y[2] != 51 || y[3] != 90 ||
            y[4] != 81 || y[5] != 76 || y[6] != 47)
            SLS_ERR("failed!");
    }
    
    // qref() and find() for diagonal elements
    // fix()
    {
        VecInt diag(7);
        diag[0] = 5; diag[1] = 1; diag[2] = 5; diag[3] = 1;
        diag[4] = 5; diag[5] = 1; diag[6] = 5;
        VecLong ind(7);
        ind[0] = 4; ind[1] = 9; ind[2] = 13; ind[3] = 18;
        ind[4] = 22; ind[5] = 27; ind[6] = 31;
        for (Long i = 0; i < 7; ++i) {
            Long ind0 = a.find(i, i);
            if (ind0 != ind[i])
                SLS_ERR("failed!");
            if (a(ind0) != diag[i])
                SLS_ERR("failed!");
            if (i % 2 == 1) {
                if (a(ind0 - 1) != 0)
                    SLS_ERR("failed!");
            }
        }

        for (Long i = 0; i < 7; ++i) {
            Long ind0 = a.find(i, i);
            a(ind0) = 100;
        }
        for (Long i = 0; i < 7; ++i) {
            if (a(i, i) != 100)
                SLS_ERR("failed!");
        }
    }

    // norm_inf()
    {
        if (norm_inf(a) != 120)
            SLS_ERR("failed!");
    }
}
