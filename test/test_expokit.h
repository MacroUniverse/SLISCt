#pragma once
#include "../SLISC/expokit/expokit.h"
#include "../SLISC/arithmetic.h"
#include "../SLISC/random.h"
#include "../SLISC/disp.h"
#include "../SLISC/sparse_arith.h"

void test_expokit()
{
	using namespace slisc;

	// === params ===========
	Long N = 40; // matrix size
	Int Nbase = 20; // # of krylov basis
	Doub t = 1; // time
	Doub tol = 0; // error tol
	// ======================

	Long i, k;
	McooComp A;
	A.resize(500);
	A.reshape(N, N);

	for (k = 0; k < 10; ++k) {
		A.trim(0);
		for (i = 0; i < N; ++i) {
			A.push(2*randDoub() - 1, i, i);
		}
		for (i = 0; i < N - 1; ++i) {
			Doub val = 2 * randDoub() - 1;
			A.push(val, i, i + 1);
			A.push(val, i + 1, i);
		}
		for (i = 0; i < N - 3; ++i) {
			Doub val = 2 * randDoub() - 1;
			A.push(val, i, i + 3);
			A.push(val, i + 3, i);
		}
		for (i = 0; i < N - 5; ++i) {
			Doub val = 2 * randDoub() - 1;
			A.push(val, i, i + 5);
			A.push(val, i + 5, i);
		}
		for (i = 0; i < N - 7; ++i) {
			Doub val = 2 * randDoub() - 1;
			A.push(val, i, i + 7);
			A.push(val, i + 7, i);
		}

		VecComp x(N), y1, y2;
		linspace(x, 1, N);

		expv<'G'>(y1, A, x, t, Nbase);
		expv<'H'>(y2, A, x, t, Nbase);

		y1 -= y2;
		if (max_abs(y1) > 1e-12) error("failed!");

		TODO: use eig.h or Eigen to calculate expA and check result!
	}
}
