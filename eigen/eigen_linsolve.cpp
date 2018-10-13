#include "eigen_linsolve.h"
using namespace Eigen;

#ifdef _MSC_VER
HouseholderQrDoub::HouseholderQrDoub(MatDoub_I &a)
{
	Map<const RMatrixXd> map_a(a.ptr(), a.nrows(), a.ncols());
	qr.compute(map_a);
}

void HouseholderQrDoub::compute(MatDoub_I &a)
{
	Map<const RMatrixXd> map_a(a.ptr(), a.nrows(), a.ncols());
	qr.compute(map_a);
}

void HouseholderQrDoub::solve(VecDoub_O &x, VecDoub_I &y)
{
	Map<const VectorXd> map_y(y.ptr(), y.size());
	Map<VectorXd> map_x(x.ptr(), x.size());
	map_x = qr.solve(map_y);
}
#endif

HouseholderQrComp::HouseholderQrComp(MatComp_I &a)
{
	Map<const RMatrixXcd> map_a(a.ptr(), a.nrows(), a.ncols());
	qr.compute(map_a);
}

void HouseholderQrComp::compute(MatComp_I &a)
{
	Map<const RMatrixXcd> map_a(a.ptr(), a.nrows(), a.ncols());
	qr.compute(map_a);
}

void HouseholderQrComp::solve(VecComp_O &x, VecComp_I &y)
{
	Map<const VectorXcd> map_y(y.ptr(), y.size());
	Map<VectorXcd> map_x(x.ptr(), x.size());
	map_x = qr.solve(map_y);
}
