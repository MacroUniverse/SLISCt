#include "../SLISC/arithmatic.h"

// test +=, -=, *=, /=
void test_self_op()
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
	// TODO: vLlong /= 2; failed!

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
void test_plus_minus_times_devide()
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