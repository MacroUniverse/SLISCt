#pragma once
#include "../SLISC/flm.h"
#include "../SLISC/scalar_arith.h"
#include "../SLISC/disp.h"

void test_flm()
{
    using namespace slisc;

    // test 2 arg constructor
    {
        Flm<Int> flm(6, 10);
        if (flm.nrows() != 7) SLS_ERR("failed!");
        if (flm.ncols() != 1) SLS_ERR("failed!");
        if (flm.Lmax() != 6) SLS_ERR("failed!");
        if (flm.Mmin() != 0) SLS_ERR("failed!");
        if (flm.Mmax() != 0) SLS_ERR("failed!");
    }

    // test 4 arg constructor
    {
        Flm<Int> flm(6, -3, 3, 10);
        if (flm.nrows() != 7) SLS_ERR("failed!");
        if (flm.ncols() != 7) SLS_ERR("failed!");
        if (flm.Lmax() != 6) SLS_ERR("failed!");
        if (flm.Mmin() != -3) SLS_ERR("failed!");
        if (flm.Mmax() != 3) SLS_ERR("failed!");
    }

    // test exist()
    {
        Flm<Int> flm(6, -3, 3, 10);
        for (Long L = 0; L <= 10; ++L) {
            for (Long M = -L*2; M <= L*2; ++M) {
                if (L <= 6 && abs(M) > 3 && flm.exist(L, M)) SLS_ERR("failed!");
                if (L <= 6 && abs(M) <= min(L,Long(3)) && !flm.exist(L, M)) SLS_ERR("failed!");
                if (L > 6 && flm.exist(L, M)) SLS_ERR("failed!");
            }
        }
    }

    // test nflm()
    {
        Flm<Int> flm0(6, 10);
        if (flm0.nf() != 7) SLS_ERR("failed!");

        Flm<Int> flm1(6, -3, 3, 10);
        if (flm1.nf() != 37) SLS_ERR("failed!");
    }

    // test get()
    {
        Flm<Int> flm(6, -3, 3, 10);
        for (Long L = 0; L <= 6; ++L) {
            for (Long M = max(-L,-Long(3)); M <= min(Long(3),L); ++M) {
                if (flm.get(L, M).size() != 10) SLS_ERR("failed!");
            }
        }
    }

    // test operator=(scalar)
    {
        Flm<Int> flm(6, -3, 3, 10);
        flm = 11;
        for (Long L = 0; L <= 6; ++L) {
            for (Long M = max(-L, -Long(3)); M <= min(Long(3), L); ++M) {
                if (flm.get(L, M) != 11) SLS_ERR("failed!");
            }
        }

        flm = -5;
        for (Long L = 0; L <= 6; ++L) {
            for (Long M = max(-L, -Long(3)); M <= min(Long(3), L); ++M) {
                if (flm.get(L, M) != -5) SLS_ERR("failed!");
            }
        }
    }

    // test operator*=()
    {
        Flm<Int> flm(6, -3, 3, 10);
        flm = 11;
        for (Long L = 0; L <= 6; ++L) {
            for (Long M = max(-L, -Long(3)); M <= min(Long(3), L); ++M) {
                if (flm.get(L, M) != 11) SLS_ERR("failed!");
            }
        }

        flm *= 2;
        for (Long L = 0; L <= 6; ++L) {
            for (Long M = max(-L, -Long(3)); M <= min(Long(3), L); ++M) {
                if (flm.get(L, M) != 22) {
                    disp(flm.get(L, M));
                    SLS_ERR("failed!");
                }
            }
        }
    }
}
