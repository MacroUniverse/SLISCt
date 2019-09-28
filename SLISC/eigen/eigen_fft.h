// FFT for any length vector
// TODO: does not need to use map<>, can use pointer for FFT input directly.

#pragma once
#include "../slisc.h"
#define EIGEN_DONT_PARALLELIZE // TODO: shouldn't appear in header file
#include "Eigen/Dense"
#include "unsupported/Eigen/FFT"

namespace slisc {

class FFT
{
private:
    Eigen::FFT<Doub> fft;
public:
    inline FFT(Bool_I scaled = false);
    typedef Eigen::Matrix<Doub, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXd;
    typedef Eigen::Matrix<Comp, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign | Eigen::RowMajor> RMatrixXcd;
    inline void fwd(VecComp_O &g, VecComp_I &f);
    inline void inv(VecComp_O &f, VecComp_I &g); // 1/N factor included.
};

inline FFT::FFT(Bool_I scaled) : fft(Eigen::default_fft_impl<Doub>(), scaled ? Eigen::FFT<Doub>::Default : Eigen::FFT<Doub>::Unscaled)
{ }

inline void FFT::fwd(VecComp_O &g, VecComp_I &f)
{
    using namespace Eigen;
    g.resize(f);
    // there is a bug in eigen, using Map<const VectorXcd> will not work for fft input, so, I have to use const_cast<>() to get a non-const pointer
    Map<VectorXcd> map_f(const_cast<Comp*>(f.ptr()), f.size());
    Map<VectorXcd> map_g(g.ptr(), g.size());
    fft.fwd(map_g, map_f);
}

inline void FFT::inv(VecComp_O &f, VecComp_I &g)
{
    using namespace Eigen;
    f.resize(g);
    // there is a bug in eigen, using Map<const VectorXcd> will not work for fft input, so, I have to use const_cast<>() to get a non-const pointer
    Map<VectorXcd> map_g(const_cast<Comp*>(g.ptr()), g.size());
    Map<VectorXcd> map_f(f.ptr(), f.size());
    fft.inv(map_f, map_g);
}

} // namespace slisc
