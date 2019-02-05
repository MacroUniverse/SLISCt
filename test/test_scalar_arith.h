#pragma once
#include "../SLISC/scalar_arith.h"
#include "../SLISC/slisc.h"

void test_scalar_arith()
{
	using namespace slisc;

	// is_equiv(s, s)	
	if (!is_equiv(3, 3)) error("failed!");
	if (!is_equiv(3.1, 3.1)) error("failed!");
	{ Doub s = 3.1; const Doub &s1 = s;
	if (!is_equiv(s, s1)) error("failed!"); }
	if (is_equiv(3., 3)) error("failed!");
	if (is_equiv(float(3.1), 3.1)) error("failed!");

	// to_num(s)
	if (!is_equiv(to_num(false), 0)) error("failed!");
	if (!is_equiv(to_num(true), 1)) error("failed!");
	if (!is_equiv(to_num(Char(-125)), -125)) error("failed!");
	if (!is_equiv(to_num(Uchar(253)), 253)) error("failed!");
	if (!is_equiv(to_num(2.3), 2.3)) error("failed!");
	if (!is_equiv(to_num(Comp(3.1,3.2)), Comp(3.1,3.2))) error("failed!");

	// operator+,-,*,/ between floatig point std::complex<> and intrinsic types
	if (!is_equiv(Comp(1.1,2.2) + 3, Comp(4.1, 2.2))) error("failed!");
	if (!is_equiv(Comp(1.1, 2.2) + Uchar(3), Comp(4.1, 2.2))) error("failed!");
	if (!is_equiv(Comp(1.1, 2.2) + float(3), Comp(4.1, 2.2))) error("failed!");
	if (!is_equiv(Comp(1.1, 2.2) + 3., Comp(4.1, 2.2))) error("failed!");
	if (!is_equiv(-2 + Comp(1., 2.2), Comp(-1, 2.2))) error("failed!");
	if (!is_equiv(Char(-2) + Comp(1., 2.2), Comp(-1, 2.2))) error("failed!");

	if (!is_equiv(Comp(1.1, 2.2) - 3, Comp(-1.9, 2.2))) error("failed!");
	if (!is_equiv(Comp(1.1, 2.2) - float(3), Comp(-1.9, 2.2))) error("failed!");
	if (!is_equiv(Comp(1.1, 2.2) - 3., Comp(-1.9, 2.2))) error("failed!");
	if (!is_equiv(1 - Comp(1.5, 2.2), Comp(-0.5, -2.2))) error("failed!");
	if (!is_equiv(true - Comp(1.5, 2.2), Comp(-0.5, -2.2))) error("failed!");
	if (!is_equiv(float(1) - Comp(1.5, 2.2), Comp(-0.5, -2.2))) error("failed!");
	if (!is_equiv(1. - Comp(1.5, 2.2), Comp(-0.5, -2.2))) error("failed!");

	if (!is_equiv(Comp(1.5, 2.5) * Uchar(3), Comp(4.5, 7.5))) error("failed!");
	if (!is_equiv(Comp(1.5, 2.5) * float(3), Comp(4.5, 7.5))) error("failed!");
	if (!is_equiv(Comp(1.5, 2.5) * 3., Comp(4.5, 7.5))) error("failed!");
	if (!is_equiv(Char(3) * Comp(1.5, 2.5), Comp(4.5, 7.5))) error("failed!");

	if (!is_equiv(Comp(4.5, 7.5) / 3, Comp(1.5, 2.5))) error("failed!");
	if (!is_equiv(Comp(1.5, 4.5) / float(3.), Comp(0.5, 1.5))) error("failed!");
	if (!is_equiv(Comp(1.5, 4.5) / 3., Comp(0.5, 1.5))) error("failed!");
	if (!is_equiv(6 / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(float(6) / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(6. / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(6. / std::complex<float>(1, 1), Comp(3, -3))) error("failed!");

	// operator+,-,*,/ between two floatig point std::complex<>	
	if (!is_equiv(std::complex<float>(1, 1) + std::complex<float>(1, 1), std::complex<float>(2, 2)))
		error("failed!");
	if (!is_equiv(std::complex<float>(1, 1) - Comp(1, 1), Comp(0, 0)))
		error("failed!");
	if (!is_equiv(std::complex<float>(1, 1) / Comp(1, 1), Comp(1, 0)))
		error("failed!");
}
