#pragma once
#include "../SLISC/expokit.h"
#include "../SLISC/arithmetic.h"
#include "../SLISC/random.h"
#include "../SLISC/disp.h"
#include "../SLISC/sparse_arith.h"
#include "../SLISC/mat_fun.h"

void test_expokit()
{
    using namespace slisc;

    // === params ===========
    Long N = 40; // matrix size
    Int Nbase = 20; // # of krylov basis
    Doub t = 1; // time
    // ======================

    Long i, k;
    McooDoub A(0,0);
    A.reserve(500);
    A.reshape(N, N);

    for (k = 0; k < 10; ++k) {
        A.trim(0);
        for (i = 0; i < N; ++i) {
            A.push(2*randDoub() - 1, i, i);
        }
        for (i = 0; i < N - 1; ++i) {
            Doub val = 2 * randDoub() - 1;
            A.push(val, i, i + 1);
            A.push(val, i + 1, i);
        }
        for (i = 0; i < N - 3; ++i) {
            Doub val = 2 * randDoub() - 1;
            A.push(val, i, i + 3);
            A.push(val, i + 3, i);
        }
        for (i = 0; i < N - 5; ++i) {
            Doub val = 2 * randDoub() - 1;
            A.push(val, i, i + 5);
            A.push(val, i + 5, i);
        }
        for (i = 0; i < N - 7; ++i) {
            Doub val = 2 * randDoub() - 1;
            A.push(val, i, i + 7);
            A.push(val, i + 7, i);
        }

        VecComp x(N), y0(N), y1(N), y2(N);
        CmatDoub A0(A.n1(), A.n2()), expA(0, 0);
        expA.resize(A0); A0 = A;
        linspace(x, 1, N);

        // reference result
        exp_mat_sym(expA, A0, t);
        mul(y0, expA, x);

        // expokit result
        y1 = x;
        expv<'G'>(y1, A, t, Nbase, norm_inf(A));
        y2 = x;
        expv<'H'>(y2, A, t, Nbase, norm_inf(A));
        
        y1 -= y0;
        if (max_abs(y1) > 5e-12)
		    SLS_ERR("failed!");
        y2 -= y0;
        if (max_abs(y2) > 5e-12)
            SLS_ERR("failed!");
    }
}
