// compile print.cpp to use in debugger
#include "../SLISC/arithmatic.h"
#include "../SLISC/time.h"

void test_print()
{
	using namespace slisc;
	VecUchar v8(3);
	linspace<Uchar>(v8, 1, 3);
	VecInt vi(3);
	linspace<Int>(vi, 1, 3);
	MatUchar A8(2, 3);
	linspace<Uchar>(A8, 1, 6);
	MatInt AI(2, 3);
	linspace<Int>(AI, 1, 6);
	Mat3Doub A3(2, 2, 2);
	linspace<Doub>(A3, 1., 8.);
	Mat3Comp C3;
	C3.resize(2, 2, 2);
	linspace<Comp>(C3, Comp(0, 1), Comp(14, 15));
	printf("set breakpoint here and try print() in debugger!\n");
	pause();
}
