// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are constants
// algorithm from Numerical Recipes 3ed
// writen by Hongyu Shi

#ifndef _TRIDAG_H_
#define _TRIDAG_H_
#include "nr3.h"

struct TridagC2
{
	Int n; // n*n matrix
	Complex lower;
	VecComplex gam, bet, diag; // U matrix data, first element not used
	TridagC2() {};
	TridagC2(VecComplex_I & diag0, const Complex upper0, const Complex lower0); // constructor
	void solve(VecComplex_O &x, VecComplex_I &b);
};

// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are real constants
// diagonal are complex constants
struct TridagC3
{
	Int n; // n*n matrix
	Complex lower;
	Complex diag;
	VecComplex gam, bet; // U matrix data, first element not used
	TridagC3() {};
	TridagC3(const Complex diag, const Complex upper, const Complex lower, const Int n0); // constructor
	void solve(VecComplex_O &x, VecComplex_I &b);
};

#endif
