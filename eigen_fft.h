// FFT for any length vector
// TODO: does not need to use map<>, can use pointer for FFT input directly.

#pragma once
#include "nr3plus.h"
#define EIGEN_DONT_PARALLELIZE
#include "Eigen/Dense"
#include "unsupported/Eigen/FFT"

class FFT
{
private:
	Eigen::FFT<Doub> fft;
public:
	FFT(Bool_I scaled = false);
	typedef Eigen::Matrix<Doub, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXd;
	typedef Eigen::Matrix<Comp, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXcd;
	void fwd(VecComp_O &g, VecComp_I &f);
	void inv(VecComp_O &f, VecComp_I &g); // 1/N factor included.
};