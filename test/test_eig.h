#include "../SLISC/eig.h"
#include "../SLISC/arithmetic.h"
#include "../SLISC/sparse_arith.h"
#include "../SLISC/random.h"
#include "../SLISC/disp.h"

void test_eig()
{
	using namespace slisc;
	{
		CmatDoub a(2, 2);
		a(0, 0) = 1.; a(1, 1) = 2.;
		a(1, 0) = 3.; a(0, 1) = 3.;
		CmatDoub eigVec(0,0); VecDoub eigVal(a.nrows()); eigVec.resize(a);
		eig_sym(eigVal, eigVec, a);

		CmatDoub eigVec1(a.nrows(), eigVec.ncols()), eigVec2(0,0);
		eigVec2.resize(eigVec);
		mul(eigVec1, a, eigVec);
		mul(eigVec2, eigVec, diag(eigVal));
		eigVec1 -= eigVec2;
		if (max_abs(eigVec1) > 1e-14) SLS_ERR("failed!");
	}
	
	CmatDoub a(50, 50);
	CmatDoub eigVec(0, 0); VecDoub eigVal(a.nrows()); eigVec.resize(a);
	CmatDoub eigVec1(0, 0), eigVec2(0, 0);
	for (Long k = 0; k < 10; ++k) {
		// fill upper triangle
		for (Long j = 0; j < 50; ++j) {
			for (Long i = 0; i < j; ++i) {
				a(i, j) = a(j, i) = 20 * rand() - 10;
			}
		}
		for (Long i = 0; i < 50; ++i)
			a(i, i) = 20 * rand() - 10;

		eig_sym(eigVal, eigVec, a);
		eigVec1.resize(a.nrows(), eigVec.ncols());
		mul(eigVec1, a, eigVec);
		eigVec2.resize(eigVec);
		mul(eigVec2, eigVec, diag(eigVal));
		eigVec1 -= eigVec2;
		auto err = max_abs(eigVec1);
		if (max_abs(eigVec1) > 5e-5) {
			cout << max_abs(eigVec1) << endl;
			SLS_ERR("failed!"); // TODO: why is windows more accurate?
		}
	}
}
