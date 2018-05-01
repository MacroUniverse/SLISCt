// test
#include <iostream>
#include "nr3.h"
using namespace std;

int main() {
	VecDoub x(3, 0.);
	MatDoub A(3, 3, 0.);
	x[1] = 1.;
	cout << x[0] << "   " << x[1] << "   " << x[2] << "   " << endl;
	A[0][1] = 1.;
	cout << A[0][0] << A[0][1] << A[0][2] << endl;
	cout << A[1][0] << A[1][1] << A[1][2] << endl;
	cout << A[2][0] << A[2][1] << A[2][2] << endl;
}

