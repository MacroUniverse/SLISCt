#pragma once

#include "SLISC/imag.h"

void test_imag()
{
	using namespace slisc;
	
	// constructors
	{
		Imag x;
		if (real(x) != 0 || x.real() != 0) error("failed!");

		Imag x1(0.25);
		if (real(x1) != 0 || x1.real() != 0) error("failed!");
		if (imag(x1) != 0.25 || x1.imag() != 0.25) error("failed!");
	}

	// operator +
	{
		// Imag + Imag
		Imag x, x1(0.25), x2(0.75);
		x = x1 + x2;
		if (imag(x) != 1) error("failed!");
	}

	// operator *
	{
		// Imag * Imag
		Doub x;
		Imag x1(0.75), x2(4);
		x = x1 * x2;
	}

	// complex<T> = imag<T>
	{
		Comp x; Imag x1(1.25);
		x = x1;
		if (real(x) != 0 || imag(x) != 1.25) error("failed!");
	}
}
