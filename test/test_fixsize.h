#pragma once
#include "../SLISC/disp.h"
#include "../SLISC/fixsize.h"
#include "../SLISC/arithmetic.h"

inline void test_fixsize()
{
	using namespace slisc;

	// === FixVec ===
	// default constructor and size()
	{
		FvecDoub<13> v;
		if (v.size() != 13) error("failed!");
		if (sizeof(v) / sizeof(Doub) != 13) error("failed!");
	}
	// constant constructor
	{
		FvecDoub<21> v(3.14);
		if (v.size() != 21) error("failed!");
		if (sizeof(v) / sizeof(Doub) != 21) error("failed!");
		if (v != 3.14) error("failed!");
	}
	// test ptr()
	{
		FvecDoub<3> v;
		Doub * p = v.ptr();
		p[0] = p[1] = p[2] = 3.1;
		if (v != 3.1) error("failed!");
	}
	// test operator[], operator(), end()
	{
		FvecDoub<3> v(0.);
		v[0] = v[1] = v[2] = 3.1;
		if (v(0) != 3.1 || v(1) != 3.1 || v(2) != 3.1) error("failed!");
		v(0) = v(1) = v(2) = 6.2;
		if (v[0] != 6.2 || v[1] != 6.2 || v[2] != 6.2) error("failed!");

		if (v.end() != 6.2 || v.end(1) != 6.2 || v.end(2) != 6.2 || v.end(3) != 6.2) error("failed!");
		v.end(1) = 2.2;
		if (v[2] != 2.2) error("failed!");
	}
	// operator=
	{
		// v = s
		FvecDoub<3> v(0.);
		v = 301;
		if (v != 301) error("failed!");

		// v = v
		FvecComp<3> v1(0.), v2(Comp(2.2, 2.2));
		v1 = v2;
		if (v1 != Comp(2.2,2.2)) error("failed!");
		v1 = v;
		if (v1 != 301.) error("failed!");
	}

	// === FixCmat ===
	// default constructor and size()
	{
		FcmatDoub<7, 9> v;
		if (v.size() != 63) error("failed!");
		if (sizeof(v) / sizeof(Doub) != 63) error("failed!");
		if (v.nrows() != 7 || v.ncols() != 9) error("failed!");
	}
	// constant constructor
	{
		FcmatDoub<3, 4> v(3.14);
		if (v.size() != 12) error("failed!");
		if (sizeof(v) / sizeof(Doub) != 12) error("failed!");
		if (v.nrows() != 3 || v.ncols() != 4) error("failed!");
		if (v != 3.14) error("failed!");
	}
	// test ptr(), operator(i), operator[j], end(), end(i)
	{
		FcmatDoub<2, 4> v(0.);
		Doub *p = v.ptr();
		p[0] = 0.; p[1] = 1.;
		v(2) = 2.; v(3) = 3.;
		v[4] = 4.; v[5] = 5.;
		v.end(2) = 6.; v.end() = 7.;
		if (v(0) != 0. || v(1) != 1. || p[2] != 2. || p[3] != 3.)
			error("failed!");
		if (p[4] != 4. || p[5] != 5.) error("failed!");
		if (p[6] != 6. || p[7] != 7.) error("failed!");
	}
	// test operator(i, j), end()
	{
		FcmatDoub<2, 2> v(0.);
		v(0, 0) = 0.; v(1, 0) = 1.;
		v(0, 1) = 2.; v(1, 1) = 3.;
		if (v[0] != 0. || v[1] != 1. || v[2] != 2. || v[3] != 3.)
			error("failed!");
	}
}
