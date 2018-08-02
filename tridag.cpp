// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are constants
// algorithm from Numerical Recipes 3ed
// writen by Hongyu Shi

#include "tridag.h"

// assuming U matrix has diagonal elements of 1
TridagC2::TridagC2(VecComp_I & diag0, const Comp upper0, const Comp lower0) :
	n{ diag0.size() }, lower{ lower0 }, diag{ diag0 }
{
	Int i;
	gam.resize(n); bet.resize(n);
	bet[0] = diag[0];
	for (i = 1; i < n; ++i) {
		gam[i] = upper0 / bet[i - 1];
		bet[i] = diag[i] - lower * gam[i];
	}
}

void TridagC2::solve(VecComp_O &x, VecComp_I &b)
{
	Long i;
	x[0] = b[0] / bet[0];
	for (i = 1; i < n; ++i) {
		x[i] = (b[i] - lower * x[i - 1]) / bet[i];
	}
	for (i = n - 2; i >= 0; --i) {
		x[i] -= gam[i + 1] * x[i + 1];
	}
}


// assuming U matrix has diagonal elements of 1
TridagC3::TridagC3(const Comp diag0, const Comp upper0, const Comp lower0, const Int n0) :
	n{ n0 }, lower{ lower0 }, diag{ diag0 }
{
	Int i;
	gam.resize(n); bet.resize(n);
	bet[0] = diag;
	for (i = 1; i < n; ++i) {
		gam[i] = upper0 / bet[i - 1];
		bet[i] = diag - lower * gam[i];
	}
}

void TridagC3::solve(VecComp_O &x, VecComp_I &b)
{
	Long i;
	x[0] = b[0] / bet[0];
	for (i = 1; i < n; ++i) {
		x[i] = (b[i] - lower * x[i - 1]) / bet[i];
	}
	for (i = n - 2; i >= 0; --i) {
		x[i] -= gam[i + 1] * x[i + 1];
	}
}