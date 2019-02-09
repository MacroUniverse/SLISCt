#include "../SLISC/arithmetic.h"

inline void test_arithmetic()
{
	using namespace slisc;

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

	// test plus(), minus(), Times(), devide()
	{
		VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
		VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
		VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

		// op(v, v, s)

		// TODO: Llong version
		vDoub1 = 1.;
		Plus(vDoub, vDoub1, 1.);
		if (vDoub != 2.) error("failed!");
		Minus(vDoub, vDoub1, 1.);
		if (vDoub != 0.) error("failed!");
		Times(vDoub, vDoub1, 2.);
		if (vDoub != 2.) error("failed!");
		Divide(vDoub, vDoub1, 2.);
		if (vDoub != 0.5) error("failed!");
		// TODO: Comp version
		vDoub = 1.;
		Plus(vComp, vDoub1, Comp(1., 1.));
		if (vComp != Comp(2., 1.)) error("failed!");

		// op(v, s, v)

		// TODO: Llong version
		vDoub1 = 1.;
		Plus(vDoub, 1., vDoub1);
		if (vDoub != 2.) error("failed!");
		Minus(vDoub, 1., vDoub1);
		if (vDoub != 0.) error("failed!");
		Times(vDoub, 2., vDoub1);
		if (vDoub != 2.) error("failed!");
		vDoub1 = 2.;
		Divide(vDoub, 2., vDoub1);
		if (vDoub != 1.) error("failed!");
		// TODO: Comp version
		vDoub1 = 1.;
		Plus(vComp, vDoub1, Comp(1., 1.));
		if (vComp != Comp(2., 1.)) error("failed!");

		// op(v, v, v)

		vLlong1 = 4; vLlong2 = 2;
		Plus(vLlong, vLlong1, vLlong2);
		if (vLlong != 6) error("failed!");
		Minus(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) error("failed!");
		Times(vLlong, vLlong1, vLlong2);
		if (vLlong != 8) error("failed!");
		Divide(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) error("failed!");

		vDoub1 = 1.; vDoub2 = 2.;
		Plus(vDoub, vDoub1, vDoub2);
		if (vDoub != 3.) error("failed!");
		Minus(vDoub, vDoub2, vDoub1);
		if (vDoub != 1.) error("failed!");
		Times(vDoub, vDoub1, vDoub2);
		if (vDoub != 2.) error("failed!");
		Divide(vDoub, vDoub1, vDoub2);
		if (vDoub != 0.5) error("failed!");

		vComp1 = Comp(1., 1.); vComp2 = Comp(2., 2.);
		Plus(vComp, vComp1, vComp2);
		if (vComp != Comp(3., 3.)) error("failed!");
		Minus(vComp, vComp2, vComp1);
		if (vComp != Comp(1., 1.)) error("failed!");
		vComp2 = 2.;
		Times(vComp, vComp1, vComp2);
		if (vComp != Comp(2., 2.)) error("failed!");
		Divide(vComp, vComp1, vComp2);
		if (vComp != Comp(0.5, 0.5)) error("failed!");
	}

	// dot product
	{
		VecComp x; VecDoub y;
		linspace(x, Comp(1.1, 1.1), Comp(3.3, 3.3), 3);
		linspace(y, 1., 3., 3);
		auto s = dot(x, y);
		if (abs(s - Comp(15.4, -15.4)) > 1e-14) error("failed!");
	}

	{
		VecDoub x; VecChar y;
		linspace(x, 1.1, 3.3, 3);
		linspace(y, 1, 3, 3);
		auto s = dot(x, y);
		if (abs(s -15.4) > 1e-14) error("failed!");
	}

	// max_abs
	{
		Vector<Fcomp> x(4); x[0] = Fcomp(3, 3); x[1] = Fcomp(5, 12); x[2] = Fcomp(-1, -1);
		if (!is_equiv(max_abs(x), 13.f)) error("failed!");
		VecComp y(4); y[0] = Comp(3, 3); y[1] = Comp(5, 12); y[2] = Comp(-1, -1);
		if (!is_equiv(max_abs(y), 13.)) error("failed!");
	}

	// matrix-vector multiplication
	{
		// MatComp 
	}
}
