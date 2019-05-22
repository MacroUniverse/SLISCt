#pragma once
#include "global.h"
#include "cmat.h"

namespace slisc {

inline void inv_mat(CmatDoub_I A)
{
#ifdef SLS_CHECK_SHAPE
	if (A.nrows() != A.ncols())
		SLS_ERR("wrong shape!");
#endif
	Long N = A.nrows();
	VecLong ipiv(N);
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, A.ptr(), N, ipiv);
	LAPACKE_dgetri(LAPACK_COL_MAJOR, N, A.ptr(), N, ipiv);
}
	
} // namespace slisc
