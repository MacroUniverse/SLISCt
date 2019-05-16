// compile print.cpp to use in debugger
#include "../SLISC/arithmetic.h"
#include "../SLISC/time.h"

void test_print()
{
	using namespace slisc;
	VecChar v8(3);
	linspace(v8, 1, 3);
	VecInt vi(3);
	linspace(vi, 1, 3);
	MatChar A8(2, 3);
	linspace(A8, 1, 6);
	MatInt AI(2, 3);
	linspace(AI, 1, 6);
	Mat3Doub A3(2, 2, 2);
	linspace(A3, 1., 8.);
	Mat3Comp C3(2, 2, 2);
	linspace(C3, Comp(0, 1), Comp(14, 15));
	printf("set breakpoint here and try print() in debugger!\n");
	pause();
}
