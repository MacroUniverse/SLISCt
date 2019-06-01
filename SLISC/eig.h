// solve eigen problem
#pragma once
#include "cmat.h"
#include <mkl.h>

namespace slisc {

// only upper triangle is needed
// eigen value in ascending order
void eig_sym(VecDoub_O eigVal, CmatDoub_O eigVec, CmatDoub_I A)
{
#ifdef SLS_CHECK_BOUNDS
	if (A.n1() != A.n2() || !shape_cmp(eigVec, A)
		|| eigVal.size() != eigVec.n1())
		SLS_ERR("wrong shape!");
#endif
	eigVec = A;
	Long N = A.n2();
	eigVal.resize(N);
	Int ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U',N, eigVec.ptr(), N, eigVal.ptr());
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
	Long N = A.n2();
	eigVal.resize(N);
	Int ret = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', N,
		(MKL_Complex16 *)eigVec.ptr(), N, eigVal.ptr());
	if (ret != 0)
		SLS_ERR("failed!");
}

} // namespace slisc
