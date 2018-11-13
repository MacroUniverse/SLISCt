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
//#include "test/test_coulomb.h"
#include "test/test_mparith.h"

using std::cout; using std::endl; //using std::conj;
void test()
{
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
	test_print();*/
	cout << "test_mparith()" << endl;
	test_mparith();
	cout << "done testing!" << endl;
}
