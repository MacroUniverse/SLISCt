// evaluate (radial) coulomb wavefunctions and their derivatives
// Note that F, dF, G, dG could be calculated simultaneously with no extra cost
// precision is currently 10 digits, set precision in "cwfcomp.h"
// This is much faster than Matlab's "hypergeom()" (about 1000-10000 times faster).

#include "cwfcomp/cwfcomp.h"

namespace slisc {
using cwfcomp::Coulomb_wave_functions;

// for scalar
inline Doub coulomb1(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
	Comp F, dF;
	cwfcomp::Coulomb_wave_functions(true, l, Z/k).F_dF(k*r, F, dF);
	return real(F);
}

// for vector/matrix and tensor
inline void coulomb10(Vbase<double> &F, Int_I l, Doub_I k, const Vbase<double> &r, Doub_I Z = -1.)
{
	Long i;
	Comp F1, dF1;
	cwfcomp::Coulomb_wave_functions f(true, l, Z/k);
	for (i = 0; i < r.size(); ++i) {
		f.F_dF(k*r(i), F1, dF1);
		F(i) = real(F1);
	}
}

inline void coulomb1(VecDoub &F, Int_I l, Doub_I k, VecDoub_I &r, Doub_I Z = -1.)
{ F.resize(r); coulomb10(F, l, k, r, Z); }

inline void coulomb1(MatDoub &F, Int_I l, Doub_I k, MatDoub_I &r, Doub_I Z = -1.)
{ F.resize(r); coulomb10(F, l, k, r, Z); }

inline void coulomb1(Mat3Doub &F, Int_I l, Doub_I k, Mat3Doub_I &r, Doub_I Z = -1.)
{ F.resize(r); coulomb10(F, l, k, r, Z); }

// for scalar
inline Doub dcoulomb1(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
	Comp F, dF;
	cwfcomp::Coulomb_wave_functions(true, l, Z/k).F_dF(k*r, F, dF);
	return real(dF);
}

// for vector/matrix and tensor
inline void dcoulomb10(Vbase<double> &dF, Int_I l, Doub_I k, const Vbase<double> &r, Doub_I Z = -1.)
{
	Long i;
	Comp F1, dF1;
	cwfcomp::Coulomb_wave_functions f(true, l, Z/k);
	for (i = 0; i < r.size(); ++i) {
		f.F_dF(k*r(i), F1, dF1);
		dF(i) = k*real(dF1);
	}
}

// derivative with respect to "r", not "k*r"
inline void dcoulomb1(VecDoub &dF, Int_I l, Doub_I k, VecDoub_I &r, Doub_I Z = -1.)
{ dF.resize(r); dcoulomb10(dF, l, k, r, Z); }

inline void dcoulomb1(MatDoub &dF, Int_I l, Doub_I k, MatDoub_I &r, Doub_I Z = -1.)
{ dF.resize(r); dcoulomb10(dF, l, k, r, Z); }

inline void dcoulomb1(Mat3Doub &dF, Int_I l, Doub_I k, Mat3Doub_I &r, Doub_I Z = -1.)
{ dF.resize(r); coulomb10(dF, l, k, r, Z); }

} // namespace slisc
