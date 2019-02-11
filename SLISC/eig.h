// solve eigen problem

#include "cmat.h"
#include <mkl.h>

namespace slisc {

// only upper triangle is needed
void eigSym(VecDoub_O eigVal, CmatDoub_O eigVec, CmatDoub_I A)
{
#ifdef SLS_CHECK_BOUNDS
	if (eigVec.nrows() != eigVec.ncols())
		error("must be a square matrix!");
#endif
	eigVec = A;
	Long N = A.ncols();
	eigVal.resize(N);
	Int ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U',N, eigVec.ptr(),N, eigVal.ptr());
	if (ret != 0) error("failed!");
}

}
