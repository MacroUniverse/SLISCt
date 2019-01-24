#include "../SLISC/input.h"

void test_input()
{
	using namespace slisc;
	using namespace std;
	Int i;
	Input inp;
	VecInt v;
	VecDoub x, y;
	Int N = 3;
	v.resize(N); v = 3;
	x.resize(N); y.resize(N);
	inp.newfile("input_test.txt");

	for (i = 0; i < N; ++i) {
		v(i) = inp.Bool("input a bool");
		inp.num2(x(i), y(i), "input 2 numbers");
	}

	cout << endl << endl;

	N += 2;
	v.resize(N); v = 3;
	x.resize(N); y.resize(N);
	inp.openfile("input_test.txt");
	for (i = 0; i < N; ++i) {
		v(i) = inp.Bool("input a bool");
		inp.num2(x(i), y(i), "input 2 numbers");
	}

	cout << endl << endl;

	v = 3;
	x.resize(N); y.resize(N);
	inp.openfile("input_test.txt");
	for (i = 0; i < N; ++i) {
		v(i) = inp.Bool("input a bool");
		inp.num2(x(i), y(i), "input 2 numbers");
	}
}
