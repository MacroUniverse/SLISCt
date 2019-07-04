// object for solving complex tridiagonal linear systems for mutiple times
// off diagonals are constants
// algorithm from Numerical Recipes 3ed
// writen by Hongyu Shi

#include "tridag.h"

// tridiag matrix * vector
void Tri2Mul(VecComp_O y, VecComp_I &diag, Comp_I upper, Comp_I lower, VecComp_I x)
{
	Int i, N = diag.size();
	if (y.size() != N) y.resize(N);
	y[0] = diag[0] * x[0] + upper * x[1];
	for (i = 0; i < N - 1; ++i)
		y[i] = lower * x[i - 1] + diag[i] * x[i] + upper * x[i + 1];
	y[N - 1] = lower * x[N - 2] + diag[N - 1] * x[N - 1];
}

// assuming U matrix has diagonal elements of 1
TridagC2::TridagC2(VecComp_I &diag0, Comp_I upper0, Comp_I lower0) :
	n{ diag0.size() }, lower{ lower0 }
{
	Int i;
	gam.resize(n); bet.resize(n);
	bet[0] = diag0[0];
	for (i = 1; i < n; ++i) {
		gam[i] = upper0 / bet[i - 1];
		bet[i] = diag0[i] - lower * gam[i];
	}
}

void TridagC2::solve(Comp *x, Comp_I *b)
{
	Int i;
	x[0] = b[0] / bet[0];
	for (i = 1; i < n; ++i) 
		x[i] = (b[i] - lower * x[i - 1]) / bet[i];
	for (i = n - 2; i >= 0; --i)
		x[i] -= gam[i + 1] * x[i + 1];
}

// b stride version
void TridagC2::solve(Comp *x, Comp_I *b, Int_I bstride)
{
	Int i, ind;
	x[0] = b[0] / bet[0];
	for (i = 1; i < n; ++i)
		x[i] = (b[bstride*i] - lower * x[i - 1]) / bet[i];
	for (i = n - 2; i >= 0; --i)
		x[i] -= gam[i + 1] * x[i + 1];
}

// assuming U matrix has diagonal elements of 1
TridagC3::TridagC3(Comp_I diag0, Comp_I upper0, Comp_I lower0, Int_I n0) :
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
	Int i;
	x[0] = b[0] / bet[0];
	for (i = 1; i < n; ++i) {
		x[i] = (b[i] - lower * x[i - 1]) / bet[i];
	}
	for (i = n - 2; i >= 0; --i) {
		x[i] -= gam[i + 1] * x[i + 1];
	}
}