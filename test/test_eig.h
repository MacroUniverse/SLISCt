#include "../SLISC/arithmetic.h"
#include "../SLISC/eig.h"
#include "../SLISC/disp.h"

void test_eig()
{
	using namespace slisc;
	CmatDoub a(2, 2);
	a(0, 0) = 1.; a(1, 1) = 2.;
	a(1, 0) = 3.; a(0, 1) = 3.;
	CmatDoub eigVec; VecDoub eigVal;
	eigSym(eigVal, eigVec, a);

	CmatDoub eigVec1, eigVec2;
	mul(eigVec1, a, eigVec);
	mul(eigVec2, eigVec, diag(eigVal));
	eigVec1 -= eigVec2;
	if (max_abs(eigVec1) > 1e-14) error("failed!");
}
