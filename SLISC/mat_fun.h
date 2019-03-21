// functions of square matrix
#pragma once
#include "arithmetic.h"
#include "sparse_arith.h"
#include "eig.h"

namespace slisc {

// out = exp(a*t)
void exp_mat_sym(CmatDoub_O out, CmatDoub_I a, Doub_I t)
{
#ifdef SLS_CHECK_SHAPE
	if (a.nrows() != a.ncols() || !shape_cmp(out, a))
		error("not a square matrix!");
#endif
	VecDoub eigVal(a.nrows()); CmatDoub eigVec; eigVec.resize(a);
	eig_sym(eigVal, eigVec, a);
	eigVal *= t;
	exp(eigVal, eigVal);
	CmatDoub temp;
	mul(temp, eigVec, diag(eigVal));
	trans(eigVec);
	mul(out, temp, eigVec);
}

}
