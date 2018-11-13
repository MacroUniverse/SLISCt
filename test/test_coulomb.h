#include "SLISC/arithmatic.h"
#include "SLISC/coulomb.h"
#include "SLISC/time.h"
#include "SLISC/disp.h"

void test_coulomb()
{
	using std::cout; using std::endl; using std::conj;
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
