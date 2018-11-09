// comprehensive test of SLISC

//#include "test/test_time.h"
//#include "test/test_slisc.h"
//#include "test/test_arithmatic.h"
//#include "test/test_interp1.h"
//#include "test/test_fft.h"
//#include "test/test_random.h"
//#include "test/test_eigen_basics.h"
//#include "test/test_eigen_linsolve.h"
//#include "test/test_eigen_fft.h"
//#include "test/test_disp.h"
//#include "test/test_print.h"

#include "SLISC/arithmatic.h"
#include "SLISC/coulomb.h"
#include "SLISC/time.h"
#include "SLISC/disp.h"

using std::cout; using std::endl; using std::conj;
void test()
{
	using namespace slisc;
	Int l = 0;
	Doub k = 2.;
	Timer time;

	Long N = 500;
	VecDoub r(2), F(2); linspace(r, 0., 10., N);
	//r(0) = 5.671342685370742;
	//r(1) = 5.691382765531062;

	time.tic();
	coulombF(F, l, k, r);
	cout << "time/eval = " << time.toc()/N << endl;
	//disp(r, 15);
	//disp(F, 15);
}

int main()
{
	// temporary test
	test();

	// systematic tests
	/*cout << "test_time()" << endl;
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
	cout << "test_disp()" << endl;
	test_disp();
	cout << "test_print()" << endl;
	test_print();
	cout << "done testing!" << endl;*/
}
