#pragma once
#include "expokit/zgexpv.h"
#include "expokit/zhexpv.h"
#include "cmatobd.h"

namespace slisc {

// matrix / vector multiplication

template <class T1, class T2, SLS_IF(is_promo<T2, T1>())>
void mul(T2 *y, const MatCoo<T1> &a, T2 *x)
{
    mul_v_coo_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

template <class T1, class T2, SLS_IF(is_promo<T2, T1>())>
void mul(T2 *y, const MatCooH<T1> &a, T2 *x)
{
    mul_v_cooh_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

template <class T1, class T2, SLS_IF(is_promo<T2, T1>())>
void mul(T2 *y, const CmatObd<T1> &a, T2 *x)
{
    mul_v_cmatobd_v(y, x, a.ptr(), a.n0(), a.nblk(), a.n1());
}

// expv()
// this function is extremely slow when used in a loop! due to dynamic memory allocation
// use ZGEXPV() for MatCoo<>, ZHEXPV() for MatCooH<>
// v cannot be empty!
template <Char Option = 0, class Tvec, class Tmat, SLS_IF(
    is_dense_vec<Tvec>() &&
    (is_Comp<contain_type<Tvec>>() || is_Doub<contain_type<Tvec>>()) &&
    (is_MatCoo<Tmat>() || is_MatCooH<Tmat>()) &&
    (is_Comp<contain_type<Tmat>>() || is_Doub<contain_type<Tmat>>())
)>
inline void expv(Tvec &v, const Tmat &mat, Doub_I t, Int_I Nkrylov, Doub_I mat_norm, Doub_I tol = 0)
{
#ifdef SLS_CHECK_SHAPE
    if (mat.n1() != mat.n2() || mat.n2() != v.size())
        SLS_ERR("wrong shape!");
#endif
    Int iflag;
    VecComp wsp(MAX(Long(10), SQR(mat.n1()*(Nkrylov + 2) + 5 * (Nkrylov + 2)) + 7));
    VecInt iwsp(MAX(Nkrylov + 2, 7));

    if (Option == 'G' || (Option == 0 && is_MatCoo<Tmat>())) {
        ZGEXPV((Int)v.size(), Nkrylov, t, v.ptr(),
            tol, mat_norm, wsp.ptr(), (Int)wsp.size(),
            iwsp.ptr(), (Int)iwsp.size(), mat, 0, iflag);
    }
    else if (Option == 'H' || (Option == 0 && is_MatCooH<Tmat>())) {
        ZHEXPV((Int)v.size(), Nkrylov, t, v.ptr(),
            tol, mat_norm, wsp.ptr(), (Int)wsp.size(),
            iwsp.ptr(), (Int)iwsp.size(), mat, 0, iflag);
    }
    else SLS_ERR("unknown!");
}

} // namespace slisc
