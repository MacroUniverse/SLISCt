#pragma once
#include "../slisc.h"
#define EIGEN_DONT_PARALLELIZE // TODO: think about how to move this to source file
#include "Eigen/Dense"
#include "Eigen/LU"

// simultaneously using HouseholderQR<MatrixXd>::solve() and HouseholderQR<RMatrixXcd>::solve() will
// cause a bug in gcc
#ifdef _MSC_VER
class HouseholderQrDoub
{
private:
	typedef Eigen::Matrix<Doub, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXd;
	Eigen::HouseholderQR<RMatrixXd> qr;
public:
	HouseholderQrDoub() = default;
	HouseholderQrDoub(MatDoub_I &a);
	void compute(MatDoub_I &a);
	void solve(VecDoub_O &x, VecDoub_I &y);
};
#endif

class HouseholderQrComp
{
private:
	typedef Eigen::Matrix<Comp, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXcd;
	Eigen::HouseholderQR<RMatrixXcd> qr;
public:
	HouseholderQrComp() = default;
	HouseholderQrComp(MatComp_I &a);
	void compute(MatComp_I &a);
	void solve(VecComp_O &x, VecComp_I &y);
};
