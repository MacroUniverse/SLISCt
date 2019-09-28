// calculate spherical harmonics Y_{lm}(\theta,\phi) and coupled spherical harmonics
#pragma once
#include "scalar_arith.h"
#include "anglib.h"
#ifdef SLS_USE_GSL
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_coupling.h>
#endif

namespace slisc {

// same definition with Wolfram Alpha
// see http://littleshi.cn/online/SphHar.html
#ifdef SLS_USE_GSL
Comp ylm(Int_I l, Int_I m, Doub_I theta, Doub_I phi)
{
    Doub A, ret0;
    Long n;
    Doub th = mod(n, theta, 2 * PI);
    Doub ph;
    if (th > PI) {
        th = 2 * PI - th;
        ph = mod(n, phi + PI, 2 * PI);
    }
    else
        ph = mod(n, phi, 2 * PI);

    if (m >= 0) {
        A = sqrt((2 * l + 1)*factorial(l - m) / (4 * PI * factorial(l + m)));
        ret0 = A * gsl_sf_legendre_Plm(l, m, cos(th));
    }
    else {
        A = pow(-1, m) * sqrt((2 * l + 1)*factorial(l + m) / (4 * PI * factorial(l - m)));
        ret0 = A * gsl_sf_legendre_Plm(l, -m, cos(th));
    }

    if (m != 0) {
        Comp ret = ret0 * exp(Comp(0., m*ph));
        return ret;
    }
    return ret0;
}
#endif

// generalized spherical harmonics
Comp yl1l2LM(Int_I l1, Int_I l2, Int_I L, Int_I M,
    Doub_I theta1, Doub_I phi1, Doub_I theta2, Doub_I phi2)
{
    Comp sum = 0.;
    Long Ndim, m1_max;
    cgTableDim(Ndim, m1_max, l1, l2, M);
    for (Long m1 = m1_max + 1 - Ndim; m1 <= m1_max; ++m1) {
        Long m2 = M - m1;
        sum += cleb(2*l1, 2*m1, 2*l2, 2*m2, 2*L, 2*M) *
            ylm(l1, m1, theta1, phi1) * ylm(l2, m2, theta2, phi2);
    }
    return sum;
}

} // namespace slisc
