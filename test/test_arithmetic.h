#include "../SLISC/arithmetic.h"

inline void test_arithmetic()
{
	using namespace slisc;

	// test scalar functions
	{
		if (!isodd(3) || !isodd(-3) || isodd(4) || isodd(-4) )
			error("failed!");
		if (!ispow2(4) || !ispow2(32) || !ispow2(1024) || !ispow2(65536)
			|| ispow2(12) || ispow2(48))
			error("failed!");
		if (mod(-1, 4) != 3 || mod(-2, 4) != 2 || mod(-3, 4) != 1 || mod(-4, 4) != 0)
			error("failed!");
	}

	// test shape_cmp
	{
		if (shape_cmp(MatComp(3,4), VecInt(12))
			|| !shape_cmp(Mat3Doub(7, 3, 5), Mat3Comp(7, 3, 5))
			|| !shape_cmp(MatDoub(3, 4), CmatDoub(3, 4))
			|| !shape_cmp(CmatInt(3, 4), MatChar(3, 4)))
			error("failed!");
	}

	// test +=, -=, *=, /=
	{
		VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
		VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
		VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

		// vector ?= scalar

		vLlong = 1;
		vLlong += 1;
		if (vLlong != 2) error("failed!");
		vLlong -= 1;
		if (vLlong != 1) error("failed!");
		vLlong *= 2;
		if (vLlong != 2) error("failed!");
		vLlong /= 2;
		if (vLlong != 1) error("failed!");

		vDoub = 1.;
		vDoub += 1.;
		if (vDoub != 2.) error("failed!");
		vDoub -= 1.;
		if (vDoub != 1.) error("failed!");
		vDoub *= 2.;
		if (vDoub != 2.) error("failed!");
		vDoub /= 2.;
		if (vDoub != 1.) error("failed!");

		vComp = Comp(1., 1.);
		vComp += Comp(1., 1.);
		if (vComp != Comp(2., 2.)) error("failed!");
		vComp -= Comp(1., 1.);
		if (vComp != Comp(1., 1.)) error("failed!");
		vComp *= 2.;
		if (vComp != Comp(2., 2.)) error("failed!");
		vComp /= 2.;
		if (vComp != Comp(1., 1.)) error("failed!");

		// vector ?= vector

		vLlong = 1; vLlong1 = 1;
		vLlong += vLlong1;
		if (vLlong != 2) error("failed!");
		vLlong -= vLlong1;
		if (vLlong != 1) error("failed!");
		vLlong2 = 2;
		vLlong *= vLlong2;
		if (vLlong != 2) error("failed!");
		vLlong /= vLlong2;
		if (vLlong != 1) error("failed!");

		vDoub = 1.; vDoub1 = 1.;
		vDoub += vDoub1;
		if (vDoub != 2.) error("failed!");
		vDoub -= vDoub1;
		if (vDoub != 1.) error("failed!");
		vDoub2 = 2.;
		vDoub *= vDoub2;
		if (vDoub != 2.) error("failed!");
		vDoub /= vDoub2;
		if (vDoub != 1.) error("failed!");

		vComp = Comp(1., 1.); vComp1 = Comp(1., 1.);
		vComp += vComp1;
		if (vComp != Comp(2., 2.)) error("failed!");
		vComp -= vComp1;
		if (vComp != Comp(1., 1.)) error("failed!");
		vComp2 = 2.;
		vComp *= vComp2;
		if (vComp != Comp(2., 2.)) error("failed!");
		vComp /= vComp2;
		if (vComp != Comp(1., 1.)) error("failed!");
	}

	// test plus(), minus(), times(), devide()
	{
		VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
		VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
		VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

		// op(v, v, s)

		// TODO: Llong version
		vDoub1 = 1.;
		plus(vDoub, vDoub1, 1.);
		if (vDoub != 2.) error("failed!");
		minus(vDoub, vDoub1, 1.);
		if (vDoub != 0.) error("failed!");
		times(vDoub, vDoub1, 2.);
		if (vDoub != 2.) error("failed!");
		divide(vDoub, vDoub1, 2.);
		if (vDoub != 0.5) error("failed!");
		// TODO: Comp version
		vDoub = 1.;
		plus(vComp, vDoub1, Comp(1., 1.));
		if (vComp != Comp(2., 1.)) error("failed!");

		// op(v, s, v)

		// TODO: Llong version
		vDoub1 = 1.;
		plus(vDoub, 1., vDoub1);
		if (vDoub != 2.) error("failed!");
		minus(vDoub, 1., vDoub1);
		if (vDoub != 0.) error("failed!");
		times(vDoub, 2., vDoub1);
		if (vDoub != 2.) error("failed!");
		vDoub1 = 2.;
		divide(vDoub, 2., vDoub1);
		if (vDoub != 1.) error("failed!");
		// TODO: Comp version
		vDoub1 = 1.;
		plus(vComp, vDoub1, Comp(1., 1.));
		if (vComp != Comp(2., 1.)) error("failed!");

		// op(v, v, v)

		vLlong1 = 4; vLlong2 = 2;
		plus(vLlong, vLlong1, vLlong2);
		if (vLlong != 6) error("failed!");
		minus(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) error("failed!");
		times(vLlong, vLlong1, vLlong2);
		if (vLlong != 8) error("failed!");
		divide(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) error("failed!");

		vDoub1 = 1.; vDoub2 = 2.;
		plus(vDoub, vDoub1, vDoub2);
		if (vDoub != 3.) error("failed!");
		minus(vDoub, vDoub2, vDoub1);
		if (vDoub != 1.) error("failed!");
		times(vDoub, vDoub1, vDoub2);
		if (vDoub != 2.) error("failed!");
		divide(vDoub, vDoub1, vDoub2);
		if (vDoub != 0.5) error("failed!");

		vComp1 = Comp(1., 1.); vComp2 = Comp(2., 2.);
		plus(vComp, vComp1, vComp2);
		if (vComp != Comp(3., 3.)) error("failed!");
		minus(vComp, vComp2, vComp1);
		if (vComp != Comp(1., 1.)) error("failed!");
		vComp2 = 2.;
		times(vComp, vComp1, vComp2);
		if (vComp != Comp(2., 2.)) error("failed!");
		divide(vComp, vComp1, vComp2);
		if (vComp != Comp(0.5, 0.5)) error("failed!");
	}

	// dot product
	{
		VecComp x; VecChar y;
		linspace(x, Comp(1.1, 1.1), Comp(3.3, 3.3), 3);
		linspace(y, 1, 3, 3);
		auto s = dot(x, y);
		if (abs(s - Comp(15.4, -15.4)) > 1e-14) error("failed!");
	}

	{
		using namespace slisc;
		using namespace std;
		Comp s1(1.1, 2.2); Comp s2 = 2;
		Comp s = s1 + s2;
		Comp s3;
	}

	{
		VecDoub x; VecChar y;
		linspace(x, 1.1, 3.3, 3);
		linspace(y, 1, 3, 3);
		auto s = dot(x, y);
		if (abs(s - Comp(15.4, -15.4)) > 1e-14) error("failed!");
	}

	// matrix-vector multiplication
	{
		// MatComp 
	}
}
