// solve eigen problem
#pragma once
#include "cmat.h"

namespace slisc {

// only upper triangle is needed
// eigen value in ascending order
template <class Tv, class Tmat,    class Tmat2, SLS_IF(
    is_dense_vec<Tv>() && is_Doub<contain_type<Tv>>() &&
    is_dense_mat<Tmat>() && is_same_major<Tmat, Tmat2>() &&
    is_dense_mat<Tmat2>() && is_Doub<contain_type<Tmat2>>())>
void eig_sym(Tv &eigVal, Tmat &eigVec, const Tmat2 &A)
{
#ifdef SLS_CHECK_BOUNDS
    if (A.n1() != A.n2() || !shape_cmp(eigVec, A)
        || eigVal.size() != eigVec.n1())
        SLS_ERR("wrong shape!");
#endif
    eigVec = A;
    Int N = (Int)A.n2();
    Int ret;
    if (is_cmajor<Tmat>())
        ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, eigVec.ptr(), N, eigVal.ptr());
    else
        ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, eigVec.ptr(), N, eigVal.ptr());
    if (ret != 0) SLS_ERR("failed!");
}

// only upper triangle is needed
// eigen value in ascending order
void eig_her(VecDoub_O eigVal, CmatComp_O eigVec, CmatComp_I A)
{
#ifdef SLS_CHECK_BOUNDS
    if (A.n1() != A.n2() || !shape_cmp(eigVec, A)
        || eigVal.size() != eigVec.n1())
        SLS_ERR("wrong shape!");
#endif
    eigVec = A;
    Int N = (Int)A.n2();
    eigVal.resize(N);
    Int ret = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', N,
        (double _Complex*)eigVec.ptr(), N, eigVal.ptr());
    if (ret != 0)
        SLS_ERR("failed!");
}

} // namespace slisc
