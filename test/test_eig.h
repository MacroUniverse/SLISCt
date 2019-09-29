#pragma once
#include "../SLISC/eig.h"
#include "../SLISC/arithmetic.h"
#include "../SLISC/sparse_arith.h"
#include "../SLISC/random.h"
#include "../SLISC/disp.h"

void test_eig()
{
    using namespace slisc;

    // test eig_sym(CmatDoub)
    {
        Long N = 20;
        CmatDoub a(N, N);
        CmatDoub eigVec(0, 0);
        VecDoub eigVal(a.n1()); eigVec.resize(a);
        CmatDoub eigVec1(0, 0), eigVec2(0, 0);
        for (Long k = 0; k < 5; ++k) {
            // fill upper triangle
            for (Long j = 0; j < N; ++j) {
                for (Long i = 0; i < j; ++i) {
                    a(i, j) = a(j, i) = 20 * rand() - 10;
                }
            }
            for (Long i = 0; i < N; ++i)
                a(i, i) = 20 * rand() - 10;

            eig_sym(eigVal, eigVec, a);
            eigVec1.resize(a.n1(), eigVec.n2());
            mul(eigVec1, a, eigVec);
            eigVec2.resize(eigVec);
            mul(eigVec2, eigVec, diag(eigVal));
            eigVec1 -= eigVec2;
            if (max_abs(eigVec1) > 5e-5) {
                cout << max_abs(eigVec1) << endl;
                SLS_ERR("failed!"); // TODO: why is Visual Studio more accurate?
            }
        }
    }
    
    // test eig_her(CmatDoub)
    {
        Long N = 20;
        CmatComp a(N, N);
        CmatComp eigVec(0, 0);
        VecDoub eigVal(a.n1()); eigVec.resize(a);
        CmatComp eigVec1(0, 0), eigVec2(0, 0);
        for (Long k = 0; k < 5; ++k) {
            // fill upper triangle
            for (Long j = 0; j < N; ++j) {
                for (Long i = 0; i < j; ++i) {
                    a(i, j) = Comp(20 * rand() - 10, 20 * rand() - 10);
                    a(j, i) = conj(a(i, j));
                }
            }
            for (Long i = 0; i < N; ++i)
                a(i, i) = 20 * rand() - 10;

            eig_her(eigVal, eigVec, a);
            eigVec1.resize(a.n1(), eigVec.n2());
            mul(eigVec1, a, eigVec);
            eigVec2.resize(eigVec);
            mul(eigVec2, eigVec, diag(eigVal));
            eigVec1 -= eigVec2;
            if (max_abs(eigVec1) > 5e-5) {
                cout << max_abs(eigVec1) << endl;
                SLS_ERR("failed!"); // TODO: why is Visual Studio more accurate?
            }
        }
    }
}
