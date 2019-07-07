// evaluate (radial) coulomb wavefunctions and their derivatives
// Note that F, dF could be calculated simultaneously with no extra cost
// TODO : implement G, dG
#pragma once
#include "slisc.h"
#include "cwfcomp/cwfcomp.h"

namespace slisc {
using cwfcomp::Coulomb_wave_functions;

// === coulombF() ===
// efficiency is about 1e-4s/evaluation
// max abs error is 4e-9.

// for scalar
inline Doub coulombF(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
	Comp F, dF;
	cwfcomp::Coulomb_wave_functions f(true, l, Z/k);
	f.F_dF(k*r, F, dF);
	return real(F);
}

// for vector/matrix and tensor
template <class Tv1, class Tv2, SLS_IF(is_dense<Tv1>() && is_dense<Tv2>())>
inline void coulombF(Tv1 &F, Int_I l, Doub_I k, const Tv2 &r, Doub_I Z = -1.)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(F, r))
		SLS_ERR("wrong shape!");
#endif

	Long i;
	Comp F1, dF1;
	Doub eta = Z / k;
	for (i = 0; i < r.size(); ++i) {
		// must use a new "Coulomb_wave_functions" to keep accuracy
		// if the same object is used, the ODE solver might solve next point based on the previous point
		// constructor time is only 1/60 to 1/10 of the calculation time, so don't try to optimize this!
		cwfcomp::Coulomb_wave_functions f(true, l, eta);
		f.F_dF(k*r(i), F1, dF1);
		F(i) = real(F1);
	}
}

// === coulombDF() ===
// accuracy not tested
// same performance

// for scalar
inline Doub coulombDF(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
	Comp F, dF;
	cwfcomp::Coulomb_wave_functions f(true, l, Z / k);
	f.F_dF(k*r, F, dF);
	return real(dF);
}

// for vector/matrix and tensor
inline void coulombDF0(Vbase<Doub> &dF, Int_I l, Doub_I k, const Vbase<Doub> &r, Doub_I Z = -1.)
{
	Long i;
	Comp F1, dF1;
	Doub eta = Z / k;
	for (i = 0; i < r.size(); ++i) {
		cwfcomp::Coulomb_wave_functions f(true, l, eta);
		f.F_dF(k*r(i), F1, dF1);
		dF(i) = real(dF1);
	}
}

// derivative with respect to "r", not "k*r"
inline void coulombDF(VecDoub &dF, Int_I l, Doub_I k, VecDoub_I r, Doub_I Z = -1.)
{ dF.resize(r); coulombDF0(dF, l, k, r, Z); }

inline void coulombDF(MatDoub &dF, Int_I l, Doub_I k, MatDoub_I r, Doub_I Z = -1.)
{ dF.resize(r); coulombDF0(dF, l, k, r, Z); }

inline void coulombDF(Mat3Doub &dF, Int_I l, Doub_I k, Mat3Doub_I r, Doub_I Z = -1.)
{ dF.resize(r); coulombDF0(dF, l, k, r, Z); }

// === coulombFDF() ===
// coulombF and derivative with no extra cost

// for scalar
inline void coulombFDF(Doub_O F, Doub_O dF, Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
	Comp F1, dF1;
	cwfcomp::Coulomb_wave_functions f(true, l, Z / k);
	f.F_dF(k*r, F1, dF1);
	F = real(F1); dF1 = real(dF1);
}

// for vector/matrix and tensor
inline void coulombFDF0(Vbase<Doub> &F, Vbase<Doub> &dF, Int_I l, Doub_I k, const Vbase<Doub> &r, Doub_I Z = -1.)
{
	Long i;
	Comp F1, dF1;
	Doub eta = Z / k;
	for (i = 0; i < r.size(); ++i) {
		// must use a new "Coulomb_wave_functions" to keep accuracy
		// if the same object is used, the ODE solver might solve next point based on the previous point
		// constructor time is only 1/60 to 1/10 of the calculation time, so don't try to optimize this!
		cwfcomp::Coulomb_wave_functions f(true, l, eta);
		f.F_dF(k*r(i), F1, dF1);
		F(i) = real(F1); dF(i) = real(dF1);
	}
}

inline void coulombFDF(VecDoub &F, Vbase<Doub> &dF, Int_I l, Doub_I k, VecDoub_I r, Doub_I Z = -1.)
{
	F.resize(r); coulombFDF0(F, dF, l, k, r, Z);
}

inline void coulombFDF(MatDoub &F, Vbase<Doub> &dF, Int_I l, Doub_I k, MatDoub_I r, Doub_I Z = -1.)
{
	F.resize(r); coulombFDF0(F, dF, l, k, r, Z);
}

inline void coulombFDF(Mat3Doub &F, Vbase<Doub> &dF, Int_I l, Doub_I k, Mat3Doub_I r, Doub_I Z = -1.)
{
	F.resize(r); coulombFDF0(F, dF, l, k, r, Z);
}

} // namespace slisc
