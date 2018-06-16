// comprehensive test of nr3plus.h

#include "nr3.h"
#include "nr3plus.h"
#include "interp_1d.h"
#include "interp_2d.h"

using namespace std;

int main()
{
	// test new disp()
	VecUchar v8(3);
	v8[0] = 1; v8[1] = 2; v8[2] = 3;
	VecInt vi(3);
	vi[0] = 1; vi[1] = 2; vi[2] = 3;
	MatUchar A8;
	A8.assign(2, 3, 0);
	A8[0][0] = 1; A8[0][1] = 3; A8[0][2] = 5; A8[1][2] = 11;
	MatInt AI;
	AI.assign(2, 3, 0);
	AI[0][0] = 1; AI[0][1] = 3; AI[0][2] = 5; AI[1][2] = 11;
	Mat3DDoub A3;
	A3.resize(2, 2, 2);
	Doub *pA3 = A3[0][0];
	for (Int i = 0; i < 8; ++i)
		pA3[i] = 1. + (Doub)i;
	Mat3DComplex C3;
	C3.resize(2, 2, 2);
	Complex *pC3 = C3[0][0];
	for (Int i = 0; i < 8; ++i)
		pC3[i] = Complex(1. + (Doub)i, (Doub)i);

	// test cubic spline interp 1D
	Int N{ 2 };
	VecDoub xx(N), yy(N);
	linspace(xx, 0., 1.); linspace(yy, 1., 2.);
	Spline_interp myfun(xx, yy);
	cout << myfun.interp(0.5) << endl;

	// test cubic spline interp 2D
	Int i, j, m{ 5 }, n{ 5 };
	MatComplex Psi(m, n), Psi0(m, n);
	VecDoub x0, y0, x, y;
	linspace(x0, 1, n, n); linspace(y0, 1, m, m);
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			Psi0[i][j] = Complex(x0[j]*x0[j], y0[n - i - 1]*y0[n - i - 1]);
	disp(Psi0);
	x = x0; y = y0;
	Spline_grid(Psi, x, y, Psi0, x0, y0);
	disp(Psi);

	VecDoub v(3,0.); VecComplex vc(3,0.);
	linspace(v, 0., PI);
	//disp(v, 8);
	//disp(v, 0, 2, 8);
	linspace(vc, 0., Complex(PI, 2.*PI));
	//disp(vc, 8);
	//disp(vc, 0, 2, 8);
	Complex c, c1(3.,5.);
	c = c1 + 1;
	cout << c << endl;
	c = c1 - 1;
	cout << c << endl;
	c = 1 - c1;
	cout << c << endl;
	c = c1 * 2;
	cout << c << endl;
	c = 2 * c1;
	cout << c << endl;
	c = 2/c1;
	cout << c << endl;
	MatDoub a(3,3,0.);
	linspace(a, 0., 8.);
	//disp(a, 5);
}
