#pragma once
#include "../SLISC/global.h"
#include "../SLISC/scalar_arith.h"

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

	// CONST<T>() CONSTC<T>()
	if (!is_equiv(Aconst<Int, 0>::value, 0)) error("failed!");
	if (!is_equiv(Aconst<Llong, 1>::value, Llong(1))) error("failed!");
	if (!is_equiv(Aconst<float, -9>::value, float(-9))) error("failed!");
	if (!is_equiv(Aconst<Doub, -1>::value, -1.)) error("failed!");
	if (!is_equiv(Aconst<Comp, 3>::value, 3.)) error("failed!");
	if (!is_equiv(Cconst<Comp, -2, -3>::value, Comp(-2,-3))) error("failed!");

	// to_num(s)
	if (!is_equiv(to_num(false), 0)) error("failed!");
	if (!is_equiv(to_num(true), 1)) error("failed!");
	if (!is_equiv(to_num(Char(-125)), -125)) error("failed!");
	if (!is_equiv(to_num(Uchar(253)), 253)) error("failed!");
	if (!is_equiv(to_num(2.3), 2.3)) error("failed!");
	if (!is_equiv(to_num(Comp(3.1,3.2)), Comp(3.1,3.2))) error("failed!");

	// test scalar functions
	if (!isodd(3) || !isodd(-3) || isodd(4) || isodd(-4))
		error("failed!");
	if (!ispow2(4) || !ispow2(32) || !ispow2(1024) || !ispow2(65536)
		|| ispow2(12) || ispow2(48))
		error("failed!");
	if (mod(-1, 4) != 3 || mod(-2, 4) != 2 || mod(-3, 4) != 1 || mod(-4, 4) != 0)
		error("failed!");

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
	if (!is_equiv(Comp(1.5, 4.5) / float(3), Comp(0.5, 1.5))) error("failed!");
	if (!is_equiv(Comp(1.5, 4.5) / 3., Comp(0.5, 1.5))) error("failed!");
	if (!is_equiv(6 / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(float(6) / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(6. / Comp(1, 1), Comp(3, -3))) error("failed!");
	if (!is_equiv(6. / std::complex<float>(1, 1), Comp(3, -3))) error("failed!");

	// operator+,-,*,/ between two floatig point std::complex<>	
	if (!is_equiv(std::complex<float>(1.25, 3.75) + std::complex<float>(0.75, -1.75), std::complex<float>(2, 2)))
		error("failed!");
	if (!is_equiv(std::complex<float>(1.25, 3.75) - Comp(-0.75, 0.75), Comp(2, 3)))
		error("failed!");
	auto z = std::complex<float>(101.5, -205.5) / Comp(101.5, -205.5);
	if (!is_equiv(std::complex<float>(1., 2.) / Comp(1., 2.), Comp(1, 0)))
		error("failed!");
}
