// partly translated from anglib.f90
// as directly as possible! all vars values unchanged.

// anglib.f90: angular momentum coupling coefficients in Fortran 90
// Copyright (C) 1998  Paul Stevenson

#pragma once
#include "scalar_arith.h"

namespace slisc {

inline Doub binom(Long_I n, Long_I r) {
    if (n == r || r == 0)
        return 1.;
    else if (r == 1)
        return (Doub)n;
    else
        return n / Doub(n - r)*binom(n - 1, r);
}

// calc dimension of CG table and max(m1)
// CGtableDim() in MyMatlabLibrary
inline void cgTableDim(Long_O Ndim, Long_O m1_max, Long_I l1, Long_I l2, Long_I M)
{
    if (abs(M) > abs(l1 - l2))
        Ndim = l1 + l2 - abs(M) + 1;
    else
        Ndim = 2 * min(l1, l2) + 1;

    if (M >= l1 - l2)
        m1_max = l1;
    else
        m1_max = l2 + M;
}

// clebsch gordan coefficient [two_j1/2,two_m1/2,two_j2/2,two_m2/2,two_j/2,two_m/2]
inline Doub cleb(Long_I two_j1, Long_I two_m1, Long_I two_j2, Long_I two_m2, Long_I two_j, Long_I two_m) {

    Long_I j1 = two_j1, m1 = two_m1, j2 = two_j2, m2 = two_m2, j = two_j, m = two_m;
    Doub cleb, factor, sum;
    Long par, z, zmin, zmax;

    if (2 * (j1 / 2) - Long(2 * (j1 / 2.0)) != 2 * (abs(m1) / 2) - Long(2 * (abs(m1) / 2.0)) ||
        2 * (j2 / 2) - Long(2 * (j2 / 2.0)) != 2 * (abs(m2) / 2) - Long(2 * (abs(m2) / 2.0)) ||
        2 * (j / 2) - Long(2 * (j / 2.0)) != 2 * (abs(m) / 2) - Long(2 * (abs(m) / 2.0)) ||
        j1<0 || j2<0  || j<0 || abs(m1)>j1 || abs(m2)>j2 ||
        abs(m)>j || j1 + j2<j || abs(j1 - j2)>j || m1 + m2 != m) {
        cleb = 0.;
    }
    else {
        factor = 0.0;
        factor = binom(j1, (j1 + j2 - j) / 2) / binom((j1 + j2 + j + 2) / 2, (j1 + j2 - j) / 2);
        factor = factor * binom(j2, (j1 + j2 - j) / 2) / binom(j1, (j1 - m1) / 2);
        factor = factor / binom(j2, (j2 - m2) / 2) / binom(j, (j - m) / 2);
        factor = sqrt(factor);

        zmin = MAX(0, max(j2 + (j1 - m1) / 2 - (j1 + j2 + j) / 2, j1 + (j2 + m2) / 2 - (j1 + j2 + j) / 2));
        zmax = min((j1 + j2 - j) / 2, min((j1 - m1) / 2, (j2 + m2) / 2));

        sum = 0.0;
        for (z = zmin; z <= zmax; ++z) {
            par = 1;
                if (2 * (z / 2) - int(2 * (z / 2.0)) != 0)
                    par = -1;
            sum = sum + par * binom((j1 + j2 - j) / 2, z)*binom((j1 - j2 + j) / 2, (j1 - m1) / 2 - z)
                * binom((-j1 + j2 + j) / 2, (j2 + m2) / 2 - z);
        }
        cleb = factor * sum;
    }
    return cleb;
}

// 3j symbol [j1/2,m1/2,j2/2,m2/2,j/2,m/2]
// written by me
inline Doub threej(Long_I j1, Long_I m1, Long_I j2, Long_I m2, Long_I j3, Long_I m3)
{
#ifndef NDEBUG
    if (isodd(j1 - j2 - m3)) SLS_ERR("unknown!");
#endif
    return pow(-1, (j1 -j2 -m3)/2) / sqrt(j3+1.) * cleb(j1, m1, j2, m2, j3, -m3);
}

inline Doub angdelta(Long_I a, Long_I b, Long_I c) {
    Doub angdelta, scr1;
    scr1 = factorial((a + b - c) / 2);
    scr1 = scr1 / factorial((a + b + c) / 2 + 1);
    scr1 = scr1 * factorial((a - b + c) / 2);
    scr1 = scr1 * factorial((-a + b + c) / 2);
    angdelta = sqrt(scr1);
    return angdelta;
}

// 6j symbol [a/2,b/2,c/2; d/2,e/2,f/2]
inline Doub sixj(Long_I a, Long_I b, Long_I c, Long_I d, Long_I e, Long_I f)
{
    Doub sixj;
    Long  nlo, nhi, n;

    Doub outfactors, sum, sumterm;

    sixj = 0.;
    if (mod(a + b, 2) != mod(c, 2)) return sixj;
    if(mod(c+d,2) != mod(e,2)) return sixj;
    if(mod(a+e,2) != mod(f,2)) return sixj;
    if(mod(b+d,2) != mod(f,2)) return sixj;
    if(abs(a-b)>c || a+b<c) return sixj;
    if(abs(c-d)>e || c+d<e) return sixj;
    if(abs(a-e)>f || a+e<f) return sixj;
    if(abs(b-d)>f || b+d<f) return sixj;

    outfactors = angdelta(a, e, f) / angdelta(a, b, c);
    outfactors = outfactors * angdelta(b, d, f)*angdelta(c, d, e);

    nlo = max( (a+b+c)/2, max((c+d+e)/2, max((b+d+f)/2, (a+e+f)/2 )));
    nhi = min( (a+b+d+e)/2, min((b+c+e+f)/2, (a+c+d+f)/2));

    sum = 0.;
    for (n = nlo; n <= nhi; ++n) {
        sumterm = pow(-1,n);
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2);
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2);
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2);
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2);
       sum = sum + sumterm;
    }

    sixj = sum * outfactors;
    return sixj;
}

// 9j symbol [a/2,b/2,c/2; d/2,e/2,f/2; g/2,h/2,i/2]
inline Doub ninej(Long_I a, Long_I b, Long_I c, Long_I d, Long_I e, Long_I f, Long_I g, Long_I h, Long_I i) {
    Doub ninej, sum;
    Long xlo, xhi, x;

    if(abs(a-b)>c || a+b<c) return 0;
    if(abs(d-e)>f || d+e<f) return 0;
    if(abs(g-h)>i || g+h<i) return 0;
    if(abs(a-d)>g || a+d<g) return 0;
    if(abs(b-e)>h || b+e<h) return 0;
    if(abs(c-f)>i || c+f<i) return 0;
    
    xlo = max(abs(b-f),max(abs(a-i),abs(h-d)));
    xhi = min(b + f, min(a + i, h + d));
    
    sum = 0.;
    for (x = xlo; x <= xhi; x += 2) {
        sum = sum + pow(-1,x)*(x + 1)*sixj(a, b, c, f, i, x)*sixj(d, e, f, b, x, h)*
            sixj(g, h, i, x, a, d);
    }
    ninej = sum;
    return ninej;
}

// calculate <y_{l1,l2}^{L,M}|y_{l,l}^{0,0}|y_{l1_,l2_}^{L_,M_}>
// implementation using 9j symbol
inline Doub yyy(Long_I l1, Long_I l2, Long_I L, Long_I M, Long_I l,
    Long_I l1_, Long_I l2_, Long_I L_, Long_I M_)
{
    Doub out = (2 * l + 1) / (4 * PI)*sqrt((2 * l1_ + 1)*(2 * l2_ + 1)*(2 * L_ + 1))*
        cleb(l * 2, 0, l1_ * 2, 0, l1 * 2, 0) *
        cleb(l * 2, 0, l2_ * 2, 0, l2 * 2, 0) *
        ninej(l * 2, l1_ * 2, l1 * 2, l * 2, l2_ * 2, l2 * 2, 0, L_ * 2, L * 2);
    return out;
}

} // namespace slisc
