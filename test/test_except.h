// test floating point exception
#include "../SLISC/scalar_arith.h"

void test_except()
{
	using namespace slisc;
	// test nan
	if (!ISNAN(NaN)) error("failed!");
	if (ISNAN(0)) error("failed!");
	
	Int i; Doub x;

	warning("comment out the following lines one by one");
	x = exp(1e10); error("failed!"); // test overflow
	x = exp(-1e10); error("failed!"); // test underflow
	x = sqrt(-1.); error("failed!"); // test invalid
	x = 1. / 0.; error("failed!"); // test division by 0
	i = 1 / 0; error("failed!"); // test integer division by 0
}
