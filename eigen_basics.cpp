#include "eigen_basics.h"

using namespace Eigen;
typedef Eigen::Matrix<Doub, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXd;
typedef Eigen::Matrix<Comp, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXcd;

void mul(MatDoub_O &v, MatDoub_I &v1, MatDoub_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (v1.ncols() != v2.nrows()) error("wrong size!")
#endif
	v.resize(v1);
	Map<const RMatrixXd> map_v1(v1.ptr(), v1.nrows(), v1.ncols()), map_v2(v2.ptr(), v2.nrows(), v2.ncols());
	Map<RMatrixXd> map_v(v.ptr(), v.nrows(), v.ncols());
	map_v = map_v1 * map_v2;
}

void mul(MatComp_O &v, MatComp_I &v1, MatComp_I &v2)
{
#ifdef _CHECKBOUNDS_
	if (v1.ncols() != v2.nrows()) error("wrong size!")
#endif
	v.resize(v1);
	Map<const RMatrixXcd> map_v1(v1.ptr(), v1.nrows(), v1.ncols()), map_v2(v2.ptr(), v2.nrows(), v2.ncols());
	Map<RMatrixXcd> map_v(v.ptr(), v.nrows(), v.ncols());
	map_v = map_v1 * map_v2;
}