// comprehensive test of SLISC

//#include "test/test_all.h"
//
//#include "SLISC/eig.h"
//#include "SLISC/arithmetic.h"
//#include "SLISC/disp.h"

//#include "test/test_meta.h"
//#include "test/test_scalar_arith.h"
//#include "test/test_arithmetic.h"
// #include "test/test_fixsize.h"
#include "test/test_sparse.h"
// using namespace slisc;

// #include <complex>
using namespace std;

//void test_eig()
//{
//	using namespace slisc;
//	CmatDoub a(2, 2);
//	a(0, 0) = 1.; a(1, 1) = 2.;
//	a(1, 0) = 3.; a(0, 1) = 3.;
//	CmatDoub eigVec; VecDoub eigVal;
//	eigSym(eigVal, eigVec, a);
//	disp(a, 10);
//	disp(eigVec, 10);
//	disp(eigVal, 10);
//	CmatDoub eigVec1, eigVec2;
//	mul(eigVec1, a, eigVec);
//	disp(eigVec1, 10);
//	mul(eigVec2, eigVec, diag(eigVal));
//	disp(eigVec2, 10);
//	eigVec1 -= eigVec2;
//	std::cout << "max error = " << max(eigVec1) << std::endl;
//}

int main()
{
	// test_all();
	/*test_meta();
	test_scalar_arith();
	test_arithmetic();*/
	//test_fixsize();
	test_sparse();
}
