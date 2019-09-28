// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are constants
// algorithm from Numerical Recipes 3ed

#pragma once

#include "nr3.h"

void Tri2Mul(VecComp_O y, VecComp_I &diag, Comp_I upper, Comp_I lower, VecComp_I x);

struct TridagC2
{
    Long n; // n*n matrix
    Comp lower;
    VecComp gam, bet; // U matrix data, first element not used
    TridagC2() {};
    TridagC2(VecComp_I & diag0, Comp_I upper0, Comp_I lower0); // constructor
    void solve(Comp *x, Comp_I *b);
    void solve(Comp *x, Comp_I *b, Int_I stride); // for matrix column, support aliasing
};

// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are real constants
// diagonal are complex constants
struct TridagC3
{
    Long n; // n*n matrix
    Comp lower;
    Comp diag;
    VecComp gam, bet; // U matrix data, first element not used
    TridagC3() {};
    TridagC3(Comp_I diag, Comp_I upper, Comp_I lower, Int_I n0); // constructor
    void solve(VecComp_O &x, VecComp_I &b);
};
