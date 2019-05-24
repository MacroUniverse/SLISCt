// solve eigen problem
#pragma once
#include "cmat.h"
#include <mkl.h>

namespace slisc {

// only upper triangle is needed
void eig_sym(VecDoub_O eigVal, CmatDoub_O eigVec, CmatDoub_I A)
{
#ifdef SLS_CHECK_BOUNDS
	if (A.nrows() != A.ncols() || !shape_cmp(eigVec, A) || eigVal.size() != eigVec.nrows())
		SLS_ERR("wrong shape!");
#endif
	eigVec = A;
	Long N = A.ncols();
	eigVal.resize(N);
	Int ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U',N, eigVec.ptr(), N, eigVal.ptr());
	if (ret != 0) SLS_ERR("failed!");
}

}
