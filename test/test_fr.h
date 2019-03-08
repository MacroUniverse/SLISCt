#pragma once
#include "../SLISC/fr.h"
#include "../SLISC/scalar_arith.h"

void test_fr()
{
	using namespace slisc;

	// test exist()
	{
		FR<Int> fr(6, -3, 3, 10);
		for (Long L = 0; L <= 10; ++L) {
			for (Long M = -L*2; M <= L*2; ++M) {
				if (L <= 6 && abs(M) > 3 && fr.exist(L, M)) error("failed!");
				if (L <= 6 && abs(M) <= min(L,Long(3)) && !fr.exist(L, M)) error("failed!");
				if (L > 6 && fr.exist(L, M)) error("failed!");
			}
		}
	}

	// test get()
	{
		FR<Int> fr(6, -3, 3, 10);
		for (Long L = 0; L <= 6; ++L) {
			for (Long M = max(-L,-Long(3)); M <= min(Long(3),L); ++M) {
				if (fr.get(L, M).size() != 10) error("failed!");
			}
		}
	}
}
