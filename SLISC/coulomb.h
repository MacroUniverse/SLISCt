// evaluate (radial) coulomb wavefunctions and their derivatives
// Note that F, dF could be calculated simultaneously with no extra cost
// TODO : implement G, dG
#pragma once
#include "global.h"
#include "cwfcomp/cwfcomp.h"
#ifdef SLS_USE_GSL
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#endif

namespace slisc {
using cwfcomp::Coulomb_wave_functions;

// === coulomb phase shift ===
#ifdef SLS_USE_GSL
Doub coulomb_sigma(Int_I l, Doub_I eta)
{
    if (eta > 0)
        SLS_ERR("are you sure?");
    gsl_sf_result lnr, arg;
    if (gsl_sf_lngamma_complex_e(l + 1, eta, &lnr, &arg) == GSL_ELOSS)
        SLS_ERR("arg has loss!");
    return arg.val;
}
#endif

// === coulombF() ===
// efficiency is about 1e-4s/evaluation
// max abs error is 4e-9.

// for scalar
inline Doub coulombF(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
    Comp F, dF;
    Doub eta = Z / k;
    cwfcomp::Coulomb_wave_functions f(true, l, eta);
    f.F_dF(k * r, F, dF);
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

    Comp F1, dF1;
    Doub eta = Z / k;
    for (Long i = 0; i < r.size(); ++i) {
        // must use a new "Coulomb_wave_functions" to keep accuracy
        // if the same object is used, the ODE solver might solve next point based on the previous point
        // constructor time is only 1/60 to 1/10 of the calculation time, so don't try to optimize this!
        cwfcomp::Coulomb_wave_functions f(true, l, eta);
        f.F_dF(k * r[i], F1, dF1);
        F[i] = real(F1);
    }
}

// === coulombDF() ===
// accuracy not tested
// same performance

// for scalar
inline Doub coulombDF(Int_I l, Doub_I k, Doub_I r, Doub_I Z = -1.)
{
    Comp F, dF;
    Doub eta = Z / k;
    cwfcomp::Coulomb_wave_functions f(true, l, eta);
    f.F_dF(k * r, F, dF);
    return real(dF);
}

// for vector/matrix and tensor
template <class Tv1, class Tv2, SLS_IF(is_dense<Tv1>() && is_dense<Tv2>())>
inline void coulombDF(Tv1 &dF, Int_I l, Doub_I k, const Tv2 &r, Doub_I Z = -1.)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(dF, r))
        SLS_ERR("wrong shape!");
#endif

    Comp F1, dF1;
    Doub eta = Z / k;
    for (Long i = 0; i < r.size(); ++i) {
        cwfcomp::Coulomb_wave_functions f(true, l, eta);
        f.F_dF(k * r[i], F1, dF1);
        dF[i] = real(dF1);
    }
}

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
