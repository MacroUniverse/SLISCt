#include "eigen_fft.h"
using Eigen::Map; using Eigen::VectorXcd;

FFT::FFT(Bool_I scaled) : fft(Eigen::default_fft_impl<Doub>(), scaled ? Eigen::FFT<Doub>::Default : Eigen::FFT<Doub>::Unscaled)
{ }

void FFT::fwd(VecComp_O &g, VecComp_I &f)
{
	g.resize(f);
	// there is a bug in eigen, using Map<const VectorXcd> will not work for fft input, so, I have to use const_cast<>() to get a non-const pointer
	Map<VectorXcd> map_f(const_cast<Comp*>(f.ptr()), f.size());
	Map<VectorXcd> map_g(g.ptr(), g.size());
	fft.fwd(map_g, map_f);
}

void FFT::inv(VecComp_O &f, VecComp_I &g)
{
	f.resize(g);
	// there is a bug in eigen, using Map<const VectorXcd> will not work for fft input, so, I have to use const_cast<>() to get a non-const pointer
	Map<VectorXcd> map_g(const_cast<Comp*>(g.ptr()), g.size());
	Map<VectorXcd> map_f(f.ptr(), f.size());
	fft.inv(map_f, map_g);
}