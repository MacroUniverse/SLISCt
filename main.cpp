// comprehensive test of SLISC

// #include "test/test_all.h"
#include "SLISC/input.h"
#include "SLISC/disp.h"

using namespace slisc;
using namespace std;

int main()
{
	// test_all();
	/*Int i;
	Input inp;
	VecInt x;
	Int N = 3;
	x.resize(N); x = 3;
	inp.newfile("test.txt");
	for (i = 0; i < N; ++i) {
		x(i) = inp.Bool("input a bool");
	}
	disp(x);

	cout << endl << endl;

	N += 2;
	x.resize(N); x = 3;
	inp.openfile("test.txt");
	for (i = 0; i < N; ++i) {
		x(i) = inp.Bool("input a bool");
	}
	disp(x);

	cout << endl << endl;

	x = 3;
	inp.openfile("test.txt");
	for (i = 0; i < N; ++i) {
		x(i) = inp.Bool("input a bool");
	}
	disp(x);*/

	Int x, y, i;
	Input in;
	in.newfile("test.txt");
	for (i = 0; i < 3; ++i)
		in.Int2_fout(x, y, "input 2 integers");
	
	TODO: generate error when input floating point number
}
