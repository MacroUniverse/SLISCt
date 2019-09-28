#pragma once
#include "../SLISC/slice_arith.h"
#include "../SLISC/disp.h"

void test_slice()
{
    using namespace slisc;
    // === test "Svector<>" ===
    CmatInt a(3,3);
    linspace(a, 1, 9);
    SvecInt v(a.ptr(), 3);
    if (v[0] != 1 || v[1] != 2 || v[2] != 3) SLS_ERR("failed!");
    v[0] = v[1] = v[2] = 0;
    if (v[0] != 0 || v[1] != 0 || v[2] != 0) SLS_ERR("failed!");
    v.set(&a(0, 1), 2);
    if (v.size() != 2) SLS_ERR("failed!");
    if (v[0] != 4 || v[1] != 5) SLS_ERR("failed!");
    v.set_ptr(&a(0, 2));
    if (v.size() != 2) SLS_ERR("failed!");
    if (v[0] != 7 || v[1] != 8 ) SLS_ERR("failed!");
    v.set_size(3);
    if (v.size() != 3) SLS_ERR("failed!");
    if (v[0] != 7 || v[1] != 8 || v[2] != 9) SLS_ERR("failed!");

    // arithmetics
    if (sum(v) != 24) SLS_ERR("failed!");
    v += 3;
    if (sum(v) !=33) SLS_ERR("failed!");

    // slice Svector<> or Dvector<> from dense container
    {
        VecInt v(10);
        linspace(v, 1, 10);
        SvecInt sv;
        slice_vec(sv, v, 2, 3);
        if (sv[0] != 3 || sv[1] != 4 || sv[2] != 5)
            SLS_ERR("failed!");
        DvecInt dv;
        slice_vec(dv, v, 2, 3, 2);
        if (dv[0] != 3 || dv[1] != 5 || dv[2] != 7)
            SLS_ERR("failed!");
    }

    // slice column from column major matrix
    {
        CmatInt a(3, 4);
        VecInt vc(3);
        SvecInt svc;
        for (Long j = 1; j < 4; ++j) {
            slice1(svc, a, j);
            vc = svc;
            for (Long i = 1; i < 3; ++i) {
                if (svc[i] != a(i, j))
                    SLS_ERR("failed!");
                if (vc[i] != svc[i])
                    SLS_ERR("failed!");
            }
        }
    }
    // slice row from row major matrix
    {
        MatInt a(3, 4); linspace(a, 1, 12);
        VecInt vr(4);
        SvecInt svr;
        for (Long i = 1; i < 3; ++i) {
            slice2(svr, a, i);
            vr = svr;
            for (Long j = 1; j < 3; ++j) {
                if (svr[j] != a(i, j))
                    SLS_ERR("failed!");
                if (vr[j] != svr[j])
                    SLS_ERR("failed!");
            }
        }
    }

    // slice row from column major matrix
    {
        CmatInt a(3, 4);
        linspace(a, 1, 12);
        VecInt v(4);
        for (Long i = 0; i < 3; ++i) {
            v = slice2(a, i);
            for (Long j = 0; j < 4; ++j) {
                if (v[j] != a(i, j))
                    SLS_ERR("failed!");
            }
            v *= 2;
            DvecInt sli; slice2(sli, a, i);
            sli = v;
            for (Long j = 0; j < 4; ++j) {
                if (v[j] != a(i, j))
                    SLS_ERR("failed!");
            }
        }

        linspace(a, 1, 12);
        DvecInt sli;
        slice2(sli, a, 0, 1, 3);
        if (sli[0] != 4 || sli[1] != 7 || sli[2] != 10)
            SLS_ERR("failed!");
        slice2(sli, a, 1, 1, 3);
        if (sli[0] != 5 || sli[1] != 8 || sli[2] != 11)
            SLS_ERR("failed!");
    }

    // slice a3(i,j,:)
    {
        Cmat3Int a3(2, 3, 4);
        linspace(a3, 1, 24);
        DvecInt sli;
        for (Long i = 0; i < 2; ++i) {
            for (Long j = 0; j < 3; ++j) {
                slice3(sli, a3, i, j);
                for (Long k = 0; k < 4; ++k) {
                    if (sli[k] != a3(i, j, k))
                        SLS_ERR("failed!");
                }
            }
        }
    }
    {
        Cmat3Int a3(2, 2, 2);
        linspace(a3, 1, 8);
        DvecInt sli;
        slice3(sli, a3, 1, 1);
        if (sli[0] != 4 || sli[1] != 8)
            SLS_ERR("failed!");
        sli /= 2;
        if (sli[0] != 2 || sli[1] != 4)
            SLS_ERR("failed!");
    }

    // slice a3(i,:,:)
    {
        Cmat3Int a3(2, 3, 4);
        linspace(a3, 1, 24);
        JcmatInt sli;
        slice23(sli, a3, 0);
        if (sli.n1() != 3 || sli.n2() != 4)
            SLS_ERR("failed!");
        if (sli.step1() != 2 || sli.step2() != 6)
            SLS_ERR("failed!");
        if (sli.size() != 12)
            SLS_ERR("failed!");
        for (Long j = 0; j < 3; ++j) {
            for (Long k = 0; k < 4; ++k) {
                if (sli(j, k) != a3(0, j, k))
                    SLS_ERR("failed!");
            }
        }
    }

    // slice a4(i,j,:,:)
    {
        Cmat4Int a4(2, 3, 4, 5);
        linspace(a4, 1, 120);
        JcmatInt sli;
        slice34(sli, a4, 1, 0);
        if (sli.n1() != 4 || sli.n2() != 5)
            SLS_ERR("failed!");
        if (sli.step1() != 6 || sli.step2() != 24)
            SLS_ERR("failed!");
        if (sli.size() != 20)
            SLS_ERR("failed!");
        for (Long k = 0; k < 4; ++k) {
            for (Long l = 0; l < 5; ++l) {
                if (sli(k, l) != a4(1, 0, k, l))
                    SLS_ERR("failed!");
            }
        }
    }

    {
        // slice Dcmat from Dcmat
        CmatInt a(6, 6);
        linspace(a, 1, 36);
        DcmatInt sa, sa2;
        slice(sa, a, 1, 4, 1, 4);
        slice(sa2, sa, 1, 2, 1, 2);
        if (sa2 != slice(a, 2, 2, 2, 2))
            SLS_ERR("failed!");

        // slice Jcmat3d from Cmat3d
        Cmat3Int a3(4, 4, 4);
        linspace(a3, 1, 64);
        Jcmat3Int s3;
        slice(s3, a3, 1, 2, 1, 2, 1, 2);
        if (s3[0] != 22 || s3[1] != 23 || s3[6] != 42 || s3[7] != 43)
            SLS_ERR("failed!");

        // slice Jcmat4d from Cmat4d
        Cmat4Int a4(4, 4, 2, 2);
        linspace(a4, 1, 64);
        Jcmat4Int s4;
        slice(s4, a4, 1, 2, 1, 2, 0, 1, 1, 1);
        if (s4[0] != 38 || s4[1] != 39 || s4[2] != 42 || s4[3] != 43)
            SLS_ERR("failed!");
    }
}
