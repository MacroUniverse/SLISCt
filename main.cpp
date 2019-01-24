// comprehensive test of SLISC

// #include "test/test_all.h"
#include "SLISC/sort.h"
#include "SLISC/random.h"
#include "SLISC/disp.h"

int main()
{
	// test_all();
	using namespace slisc;
	Long N = 10;
	VecDoub x(N);
	VecInt order(N);
	rand(x);
	linspace(order, 0, N-1);
	disp(x);
	disp(order);
	sort2(x, order);
	disp(x);
	disp(order);
}
