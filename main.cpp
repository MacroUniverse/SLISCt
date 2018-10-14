// comprehensive test of nr3plus.h

#include "test/test_time.h"
#include "test/test_slisc.h"
#include "test/test_arithmatic.h"
#include "test/test_interp1.h"
#include "test/test_fft.h"
#include "test/test_random.h"
#include "test/test_eigen_basics.h"
#include "test/test_eigen_linsolve.h"
#include "test/test_eigen_fft.h"

// === global variables ===
#include "SLISC/global.inl"

using std::cout; using std::endl; using std::conj;

// new test scratch

void test()
{
}

int main()
{
	// temporary test
	test();

	//bench_read_write();

	// systematic tests
	cout << "test_time()" << endl;
	test_time();
	cout << "test_slisc()" << endl;
	test_slisc();
	cout << "test_self_op()" << endl;
	test_self_op();
	cout << "test_plus_minus_times_devide()" << endl;
	test_plus_minus_times_devide();
	cout << "test_interp1()" << endl;
	test_interp1();
	cout << "test_fft()" << endl;
	test_fft();
	cout << "test_rand()" << endl;
	test_random();
	cout << "test_eigen_basics()" << endl;
	test_eigen_basics();
	cout << "test_eigen_linsolve()" << endl;
	test_eigen_linsolve();
	cout << "test_eigen_fft()" << endl;
	test_eigen_fft();
	cout << "done testing!" << endl;

	//// test new disp()
	//VecUchar v8(3);
	//linspace<Uchar>(v8, 1, 3);
	//VecInt vi(3);
	//linspace<Int>(vi, 1, 3);
	//MatUchar A8(2, 3);
	//linspace<Uchar>(A8, 1, 6);
	//MatInt AI(2, 3);
	//linspace<Int>(AI, 1, 6);
	//Mat3Doub A3(2, 2, 2);
	//linspace<Doub>(A3, 1., 8.);
	//Mat3Comp C3;
	//C3.resize(2, 2, 2);
	//linspace<Comp>(C3, Comp(0, 1), Comp(14, 15));
	//
	//// test cubic spline interp 1D
	//Int N{ 2 };
	//VecDoub xx(N), yy(N);
	//linspace(xx, 0., 1.); linspace(yy, 1., 2.);
	//Spline_interp myfun(xx, yy);
	//cout << myfun.interp(0.5) << endl;

	//// test cubic spline interp 2D
	//Int i, j, m{ 5 }, n{ 5 };
	//MatComp Psi(m, n), Psi0(m, n);
	//VecDoub x0, y0, x, y;
	//linspace<Doub>(x0, 1, n, n); linspace<Doub>(y0, 1, m, m);
	//for (i = 0; i < m; ++i)
	//	for (j = 0; j < n; ++j)
	//		Psi0[i][j] = Comp(x0[j]*x0[j], y0[n - i - 1]*y0[n - i - 1]);
	//disp(Psi0);
	//x = x0; y = y0;
	//Spline_grid(Psi, x, y, Psi0, x0, y0);
	//disp(Psi);

	//VecDoub v(3,0.); VecComp vc(3,0.);
	//linspace(v, 0., PI);
	////disp(v, 8);
	////disp(v, 0, 2, 8);
	//linspace<Comp>(vc, 0., Comp(PI, 2.*PI));
	////disp(vc, 8);
	////disp(vc, 0, 2, 8);
	//Comp c, c1(3.,5.);
	//c = c1 + 1;
	//cout << c << endl;
	//c = c1 - 1;
	//cout << c << endl;
	//c = 1 - c1;
	//cout << c << endl;
	//c = c1 * 2;
	//cout << c << endl;
	//c = 2 * c1;
	//cout << c << endl;
	//c = 2/c1;
	//cout << c << endl;
	//MatDoub a(3,3,0.);
	//linspace(a, 0., 8.);
	////disp(a, 5);
}
