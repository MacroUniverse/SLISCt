#include "eigen_linsolve.h"
using namespace Eigen;

typedef Matrix<Doub, Dynamic, Dynamic, AutoAlign | RowMajor> RMatrixXd;
typedef Matrix<Comp, Dynamic, Dynamic, AutoAlign | RowMajor> RMatrixXcd;

void test_eigen_linsolve()
{
	/*using std::cout;
	RMatrixXcd m(2, 2);
	m << 1., 2., 2., 4.;
	VectorXcd b(2);
	b << -1., -2.;
	VectorXcd x;

	HouseholderQR<RMatrixXcd> qr(m);
	x = qr.solve(b);
	cout << x;*/

	MatComp a(2, 2); a(0) = 1.; a(1) = 2.; a(2) = 2.; a(3) = 4.;
	VecComp y(2); y[0] = -1.; y[1] = -2.;
	VecComp x(2);
	Map<RMatrixXcd> map_a(a.ptr(), a.nrows(), a.ncols());
	Map<VectorXcd> map_y(y.ptr(), y.size());
	Map<VectorXcd> map_x(x.ptr(), x.size());

	HouseholderQR<RMatrixXcd> qr(map_a);
	map_x = qr.solve(map_y);
	
	disp(x);
}