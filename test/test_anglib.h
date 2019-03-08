#pragma once
#include "../SLISC/anglib.h"

void test_anglib()
{
	using namespace slisc;
	// test_all();
	if (factorial(5) != 120) error("failed!");
	if (factorial(6) != 720) error("failed!");
	if (binom(5, 3) != 10) error("failed!");
	if (binom(10, 6) != 210) error("failed!");

	if (abs(cleb(10, 0, 8, 0, 2, 0) - sqrt(5. / 33.)) > 1e-15) error("failed!");
	if (abs(threej(12, 0, 8, 0, 4, 0) - sqrt(5. / 143.)) > 1e-15) error("failed!");
	if (abs(sixj(2, 4, 6, 4, 2, 4) - 1. / (5 * sqrt(21.))) > 1e-15) error("failed!");
}