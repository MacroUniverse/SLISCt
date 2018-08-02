// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are constants
// algorithm from Numerical Recipes 3ed
// writen by Hongyu Shi

#ifndef _TRIDAG_H_
#define _TRIDAG_H_
#include "nr3.h"

struct TridagC2
{
	Long n; // n*n matrix
	Comp lower;
	VecComp gam, bet, diag; // U matrix data, first element not used
	TridagC2() {};
	TridagC2(VecComp_I & diag0, const Comp upper0, const Comp lower0); // constructor
	void solve(VecComp_O &x, VecComp_I &b);
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

#endif
