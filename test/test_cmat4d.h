#pragma once
#include "../SLISC/cmat4d.h"
#include "../SLISC/slice_arith.h"

void test_cmat4d()
{
	using namespace slisc;
	Cmat4Int a(2, 3, 4, 5);
	linspace(a, 1, a.size());
	if (a.dim1() != 2 || a.dim2() != 3 ||
		a.dim3() != 4 || a.dim4() != 5)
		SLS_ERR("failed!");

	Long ind = 0;
	for (Long l = 0; l < 5; ++l)
		for (Long k = 0; k < 4; ++k)
			for (Long j = 0; j < 3; ++j)
				for (Long i = 0; i < 2; ++i) {
					++ind;
					if (a(i, j, k, l) != ind)
						SLS_ERR("failed!");
				}

	for (Long l = 0; l < 5; ++l) {
		for (Long k = 0; k < 4; ++k) {
			if (slice_mat12(a, k, l).ptr() != &a(0, 0, k, l))
				SLS_ERR("failed!");
		}
	}
}
