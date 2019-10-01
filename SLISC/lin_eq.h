#pragma once
#include "global.h"
#include "cmat.h"

namespace slisc {

#ifdef SLS_USE_LAPACKE
template<class Tmat, SLS_IF(
    is_dense_mat<Tmat>() && is_Doub<contain_type<Tmat>>())>
inline void inv_mat(Tmat &A)
{
#ifdef SLS_CHECK_SHAPE
    if (A.n1() != A.n2())
        SLS_ERR("wrong shape!");
#endif
    Int N = (Int)A.n1();
    VecInt ipiv(N);
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, A.ptr(), N, ipiv.ptr());
    LAPACKE_dgetri(LAPACK_COL_MAJOR, N, A.ptr(), N, ipiv.ptr());
}
#endif
    
} // namespace slisc
