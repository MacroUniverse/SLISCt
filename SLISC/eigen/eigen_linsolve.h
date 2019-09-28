#pragma once
#include "../slisc.h"
#define EIGEN_DONT_PARALLELIZE // TODO: think about how to move this to source file
#include "Eigen/Dense"
#include "Eigen/LU"

// simultaneously using HouseholderQR<MatrixXd>::solve() and HouseholderQR<RMatrixXcd>::solve() will
// cause a bug in gcc
// TODO: now even using one will cause a bug!
#ifdef _MSC_VER

namespace slisc {

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

HouseholderQrDoub::HouseholderQrDoub(MatDoub_I &a)
{
    Eigen::Map<const RMatrixXd> map_a(a.ptr(), a.n1(), a.n2());
    qr.compute(map_a);
}

void HouseholderQrDoub::compute(MatDoub_I &a)
{
    Eigen::Map<const RMatrixXd> map_a(a.ptr(), a.n1(), a.n2());
    qr.compute(map_a);
}

void HouseholderQrDoub::solve(VecDoub_O &x, VecDoub_I &y)
{
    using namespace Eigen;
    Map<const VectorXd> map_y(y.ptr(), y.size());
    Map<VectorXd> map_x(x.ptr(), x.size());
    map_x = qr.solve(map_y);
}

HouseholderQrComp::HouseholderQrComp(MatComp_I &a)
{
    Eigen::Map<const RMatrixXcd> map_a(a.ptr(), a.n1(), a.n2());
    qr.compute(map_a);
}

void HouseholderQrComp::compute(MatComp_I &a)
{
    Eigen::Map<const RMatrixXcd> map_a(a.ptr(), a.n1(), a.n2());
    qr.compute(map_a);
}

void HouseholderQrComp::solve(VecComp_O &x, VecComp_I &y)
{
    using namespace Eigen;
    Map<const VectorXcd> map_y(y.ptr(), y.size());
    Map<VectorXcd> map_x(x.ptr(), x.size());
    map_x = qr.solve(map_y);
}

} // namespace slisc

#endif
