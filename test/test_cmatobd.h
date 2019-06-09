#pragma once
#include "../SLISC/cmatobd.h"

void test_cmatobd()
{
	using namespace slisc;

	// construct from Cmat3d
	{
		Cmat3Int a0(3, 3, 4);
		for (Long k = 0; k < 4; ++k) {
			ScmatInt sli; slice12(sli, a0, k);
			linspace(sli, 1, 9); sli.end() = 1;
		}
		CmatObd<Int> a(a0);
		VecInt x(9), y(9); linspace(x, 1, 9);
		mul(y, a, x);
		if (y[0] != 30 || y[1] != 36 || y[2] != 72 || y[3] != 66 ||
			y[4] != 116 || y[5] != 96 || y[6] != 160 || y[7] != 126 ||
			y[8] != 78)
			SLS_ERR("failed!");
	}
	
	// construct from Mcoo
	{
		Cmat3Int a0(3, 3, 4);
		for (Long k = 0; k < 4; ++k) {
			ScmatInt sli; slice12(sli, a0, k);
			linspace(sli, 1, 9); sli.end() = 1;
		}

		CmatObd<Int> a(a0);
		McooInt a1(9, 9, 36);
		for (Long i = 0; i < 9; ++i) {
			for (Long j = 0; j < 9; ++j) {
				if (a(i, j) != 0)
					a1.push(a(i, j), i, j);
			}
		}

		CmatObd<Int> a2(a1, 3, 4);

		VecInt x(9), y(9); linspace(x, 1, 9);
		mul(y, a2, x);
		if (y[0] != 30 || y[1] != 36 || y[2] != 72 || y[3] != 66 ||
			y[4] != 116 || y[5] != 96 || y[6] != 160 || y[7] != 126 ||
			y[8] != 78)
			SLS_ERR("failed!");
	}
}
