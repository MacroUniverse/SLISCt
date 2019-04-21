#include "../SLISC/arithmetic.h"
#include "../SLISC/disp.h"

inline void test_arithmetic()
{
	using namespace slisc;

	// test shape_cmp
	{
		if (shape_cmp(MatComp(3,4), VecInt(12))) SLS_ERR("failed!");
		if (!shape_cmp(Mat3Doub(7, 3, 5), Mat3Comp(7, 3, 5))) SLS_ERR("failed!");
		if (!shape_cmp(MatDoub(3, 4), CmatDoub(3, 4))) SLS_ERR("failed!");
		if (!shape_cmp(CmatInt(3, 4), MatChar(3, 4))) SLS_ERR("failed!");
		if (!shape_cmp(CmatInt(3, 4), FixCmat<Char, 3, 4>())) SLS_ERR("failed!");
	}

	// sum, max, max_abs, norm2
	{
		Long ind;
		VecBool a(4); a[0] = 1; a[1] = 0; a[2] = 1; a[3] = 1;
		if (!is_equiv(sum(a), Long(3))) SLS_ERR("failed!");

		VecInt b(4); b[0] = 2; b[1] = 3; b[2] = -1; b[3] = 5;
		if (!is_equiv(sum(b), Long(9))) SLS_ERR("failed!");
		if (!is_equiv(max(b), 5)) SLS_ERR("failed!");
		max(ind, b); if (ind != 3) SLS_ERR("failed!");
		if (!is_equiv(max_abs(b), 5)) SLS_ERR("failed!");

		VecDoub c(10); linspace(c, 1., 10.);
		if (!is_equiv(sum(c), 55.)) SLS_ERR("failed!");
		if (!is_equiv(max(c), 10.)) SLS_ERR("failed!");
		max(ind, c); if (ind != 9) SLS_ERR("failed!");
		if (!is_equiv(max_abs(c), 10.)) SLS_ERR("failed!");
		if (!is_equiv(norm2(c), 385.)) SLS_ERR("failed!");

		MatComp d(3, 3); linspace(d, Comp(1., -1.), Comp(9., -9.));
		if (!is_equiv(sum(d), Comp(45.,-45.))) SLS_ERR("failed!");
		if (!is_equiv(max_abs(d), abs(Comp(9,9)))) SLS_ERR("failed!");
		if (!is_equiv(norm2(d), 285.*2)) SLS_ERR("failed!");
	}

	// flip
	{
		VecInt v(6), v1(6), v2(6);
		linspace(v, -2, 3); v2 = v;
		linspace(v1, 3, -2);
		flip(v);
		if (v != v1) SLS_ERR("failed!");
		flip(v, v1);
		if (v != v2) SLS_ERR("failed!");
	}

	// trans
	{
		CmatComp a(3,3); linspace(a, Comp(0, 0), Comp(8, 8));
		MatComp b(3,3); linspace(b, Comp(0, 0), Comp(8, 8));
		trans(a);
		if (a != b) SLS_ERR("failed!");

		a.resize(2, 3);
		linspace(a, Comp(0, 0), Comp(5, 5));
		b.resize(3, 2);
		linspace(b, Comp(0, 0), Comp(5, 5));
		MatComp c(a.ncols(), a.nrows()); trans(c, a);
		if (c != b)  SLS_ERR("failed!");
	}

	// her
	{
		CmatComp a(3,3); linspace(a, Comp(0, 0), Comp(8, 8));
		MatComp b(3,3); linspace(b, Comp(0, 0), Comp(8, -8));
		her(a);
		if (a != b) SLS_ERR("failed!");

		a.resize(2, 3);
		linspace(a, Comp(0, 0), Comp(5, 5));
		b.resize(3, 2);
		linspace(b, Comp(0, 0), Comp(5, -5));
		MatComp c(a.ncols(), a.nrows()); her(c, a);
		if (c != b)  SLS_ERR("failed!");
	}

	// +=, -=, *=, /=
	{
		VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
		VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
		VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

		// v ?= s
		vLlong = 1;
		vLlong += 1;
		if (vLlong != 2) SLS_ERR("failed!");
		vLlong -= 1;
		if (vLlong != 1) SLS_ERR("failed!");
		vLlong *= 2;
		if (vLlong != 2) SLS_ERR("failed!");
		vLlong /= 2;
		if (vLlong != 1) SLS_ERR("failed!");

		vDoub = 1.;
		vDoub += 1.;
		if (vDoub != 2.) SLS_ERR("failed!");
		vDoub -= 1.;
		if (vDoub != 1.) SLS_ERR("failed!");
		vDoub *= 2.;
		if (vDoub != 2.) SLS_ERR("failed!");
		vDoub /= 2.;
		if (vDoub != 1.) SLS_ERR("failed!");

		vComp = Comp(1., 1.);
		vComp += Comp(1., 1.);
		if (vComp != Comp(2., 2.)) SLS_ERR("failed!");
		vComp -= Comp(1., 1.);
		if (vComp != Comp(1., 1.)) SLS_ERR("failed!");
		vComp *= 2.;
		if (vComp != Comp(2., 2.)) SLS_ERR("failed!");
		vComp /= 2.;
		if (vComp != Comp(1., 1.)) SLS_ERR("failed!");

		// v ?= v

		vLlong = 1; vLlong1 = 1;
		vLlong += vLlong1;
		if (vLlong != 2) SLS_ERR("failed!");
		vLlong -= vLlong1;
		if (vLlong != 1) SLS_ERR("failed!");
		vLlong2 = 2;
		vLlong *= vLlong2;
		if (vLlong != 2) SLS_ERR("failed!");
		vLlong /= vLlong2;
		if (vLlong != 1) SLS_ERR("failed!");

		vDoub = 1.; vDoub1 = 1.;
		vDoub += vDoub1;
		if (vDoub != 2.) SLS_ERR("failed!");
		vDoub -= vDoub1;
		if (vDoub != 1.) SLS_ERR("failed!");
		vDoub2 = 2.;
		vDoub *= vDoub2;
		if (vDoub != 2.) SLS_ERR("failed!");
		vDoub /= vDoub2;
		if (vDoub != 1.) SLS_ERR("failed!");

		vComp = Comp(1., 1.); vComp1 = Comp(1., 1.);
		vComp += vComp1;
		if (vComp != Comp(2., 2.)) SLS_ERR("failed!");
		vComp -= vComp1;
		if (vComp != Comp(1., 1.)) SLS_ERR("failed!");
		vComp2 = 2.;
		vComp *= vComp2;
		if (vComp != Comp(2., 2.)) SLS_ERR("failed!");
		vComp /= vComp2;
		if (vComp != Comp(1., 1.)) SLS_ERR("failed!");
	}

	// rem(), mod()
	{
		VecInt v(40), v1(40);
		linspace(v, -20, 19); linspace(v1, 0, 39);
		mod(v, 5); rem(v1, 5);
		if (v != v1) SLS_ERR("failed!");
	}

	{
		VecLlong v(40), v1(40);
		linspace(v, -20, 19); linspace(v1, 0, 39);
		mod(v, 5LL); rem(v1, 5LL);
		if (v != v1) SLS_ERR("failed!");
	}

	// Plus(), Minus(), Times(), Devide()
	{
		VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
		VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
		VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

		// v = v ? s
		vLlong1 = 1;
		Plus(vLlong, vLlong1, 1);
		if (vLlong != 2) SLS_ERR("failed!");
		Minus(vLlong, vLlong1, 1);
		if (vLlong != 0) SLS_ERR("failed!");
		Times(vLlong, vLlong1, 2);
		if (vLlong != 2) SLS_ERR("failed!");
		Divide(vLlong, vLlong1, 2);
		if (vLlong != 0) SLS_ERR("failed!");

		vDoub1 = 1.;
		Plus(vDoub, vDoub1, 1.);
		if (vDoub != 2.) SLS_ERR("failed!");
		Minus(vDoub, vDoub1, 1.);
		if (vDoub != 0.) SLS_ERR("failed!");
		Times(vDoub, vDoub1, 2.);
		if (vDoub != 2.) SLS_ERR("failed!");
		Divide(vDoub, vDoub1, 2.);
		if (vDoub != 0.5) SLS_ERR("failed!");
		
		vComp1 = 1.;
		Plus(vComp, vComp1, 1.);
		if (vComp != 2.) SLS_ERR("failed!");
		Minus(vComp, vComp1, 1.);
		if (vComp != 0.) SLS_ERR("failed!");
		Times(vComp, vComp1, 2.);
		if (vComp != 2.) SLS_ERR("failed!");
		Divide(vComp, vComp1, 2.);
		if (vComp != 0.5) SLS_ERR("failed!");

		// v = s ? v

		vLlong1 = 1;
		Plus(vLlong, 1, vLlong1);
		if (vLlong != 2) SLS_ERR("failed!");
		Minus(vLlong, 1, vLlong1);
		if (vLlong != 0) SLS_ERR("failed!");
		Times(vLlong, 2, vLlong1);
		if (vLlong != 2) SLS_ERR("failed!");
		vLlong1 = 2;
		Divide(vLlong, 2, vLlong1);
		if (vLlong != 1) SLS_ERR("failed!");

		vDoub1 = 1.;
		Plus(vDoub, 1., vDoub1);
		if (vDoub != 2.) SLS_ERR("failed!");
		Minus(vDoub, 1., vDoub1);
		if (vDoub != 0.) SLS_ERR("failed!");
		Times(vDoub, 2., vDoub1);
		if (vDoub != 2.) SLS_ERR("failed!");
		vDoub1 = 2.;
		Divide(vDoub, 2., vDoub1);
		if (vDoub != 1.) SLS_ERR("failed!");

		vComp1 = Comp(1.,1.);
		Plus(vComp, Comp(1., 1.), vComp1);
		if (vComp != Comp(2.,2.)) SLS_ERR("failed!");
		Minus(vComp, Comp(1.,1.), vComp1);
		if (vComp != 0.) SLS_ERR("failed!");
		Times(vComp, 2., vComp1);
		if (vComp != Comp(2.,2.)) SLS_ERR("failed!");
		vComp1 = 2.;
		Divide(vComp, Comp(2.,2.), vComp1);
		if (vComp != Comp(1.,1.)) SLS_ERR("failed!");

		// v = v ? v

		vLlong1 = 4; vLlong2 = 2;
		Plus(vLlong, vLlong1, vLlong2);
		if (vLlong != 6) SLS_ERR("failed!");
		Minus(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) SLS_ERR("failed!");
		Times(vLlong, vLlong1, vLlong2);
		if (vLlong != 8) SLS_ERR("failed!");
		Divide(vLlong, vLlong1, vLlong2);
		if (vLlong != 2) SLS_ERR("failed!");

		vDoub1 = 1.; vDoub2 = 2.;
		Plus(vDoub, vDoub1, vDoub2);
		if (vDoub != 3.) SLS_ERR("failed!");
		Minus(vDoub, vDoub2, vDoub1);
		if (vDoub != 1.) SLS_ERR("failed!");
		Times(vDoub, vDoub1, vDoub2);
		if (vDoub != 2.) SLS_ERR("failed!");
		Divide(vDoub, vDoub1, vDoub2);
		if (vDoub != 0.5) SLS_ERR("failed!");

		vComp1 = Comp(1., 1.); vComp2 = Comp(2., 2.);
		Plus(vComp, vComp1, vComp2);
		if (vComp != Comp(3., 3.)) SLS_ERR("failed!");
		Minus(vComp, vComp2, vComp1);
		if (vComp != Comp(1., 1.)) SLS_ERR("failed!");
		vComp2 = 2.;
		Times(vComp, vComp1, vComp2);
		if (vComp != Comp(2., 2.)) SLS_ERR("failed!");
		Divide(vComp, vComp1, vComp2);
		if (vComp != Comp(0.5, 0.5)) SLS_ERR("failed!");
	}

	// real(v), v = real(v)
	{
		VecComp v(3); linspace(v, Comp(1.1,1.1), Comp(3.3,3.3));
		VecComp v1(3); linspace(v1, 1.1, 3.3);
		real(v);
		if (v != v1) SLS_ERR("failed!");
		linspace(v, Comp(1.1, 1.1), Comp(3.3, 3.3));
		VecDoub v2(v.size());
		real(v2, v);
		if(v2 != v1) SLS_ERR("failed!");
	}

	// imag(v), v = imag(v)
	{
		VecComp v(3); linspace(v, Comp(1.1, 1.1), Comp(3.3, 3.3));
		VecComp v1(3); linspace(v1, 1.1, 3.3);
		imag(v);
		if (v != v1) SLS_ERR("failed!");
		linspace(v, Comp(1.1, 1.1), Comp(3.3, 3.3));
		VecDoub v2(v.size());
		imag(v2, v);
		if (v2 != v1) SLS_ERR("failed!");
	}

	// abs(v), v = abs(v)
	{
		// doub
		{
			VecDoub v(5); linspace(v, -2, 2);
			VecComp v1(5); v1[0] = 2; v1[1] = 1; v1[2] = 0; v1[3] = 1; v1[4] = 2;
			abs(v);
			if (v != v1) SLS_ERR("failed!");
			linspace(v, -2, 2);
			VecDoub v2(v1.size());
			abs(v2, v);
			if (v2 != v1) SLS_ERR("failed!");
		}
		
		// comp
		{
			VecComp v(3); linspace(v, Comp(3, 4), Comp(9, 12));
			VecComp v1(3); linspace(v1, 5, 15);
			abs(v);
			if (v != v1) SLS_ERR("failed!");
			linspace(v, Comp(3, 4), Comp(9, 12));
			VecDoub v2(v.size());
			abs(v2, v);
			if (v2 != v1) SLS_ERR("failed!");
		}
	}

	// conj(v), v = conj(v)
	{
		VecComp v(3); linspace(v, Comp(1.1, 1.1), Comp(3.3, 3.3));
		VecComp v1(3); linspace(v1, Comp(1.1, -1.1), Comp(3.3, -3.3));
		conj(v);
		if (v != v1) SLS_ERR("failed!");
		linspace(v, Comp(1.1, 1.1), Comp(3.3, 3.3));
		Vector<Lcomp> v2(v.size());
		conj(v2, v);
		if (v2 != v1) SLS_ERR("failed!");
	}

	// s = dot(v, v)
	{
		{
			VecComp x(3); VecDoub y(3);
			linspace(x, Comp(1.1, 1.1), Comp(3.3, 3.3));
			linspace(y, 1., 3.);
			auto s = dot(x, y);
			if (abs(s - Comp(15.4, -15.4)) > 1e-14) SLS_ERR("failed!");
		}

		{
			VecDoub x(3); VecChar y(3);
			linspace(x, 1.1, 3.3);
			linspace(y, 1, 3);
			auto s = dot(x, y);
			if (abs(s - 15.4) > 1e-14) SLS_ERR("failed!");
		}
	}

	// matrix-vector multiplication
	{
		CmatComp a(4,7); linspace(a, Comp(1, -1), Comp(28, -28));
		VecChar v(7); linspace(v, 1, 7);
		Vector<Lcomp> v1(4), v2(a.nrows());
		v1[0] = Comp(476, -476); v1[1] = Comp(504, -504);
		v1[2] = Comp(532, -532); v1[3] = Comp(560, -560);
		mul(v2, a, v);
		if (v2 != v1) SLS_ERR("failed!");
	}

	// vector-matrix multiplication
	{
		MatComp a(7,4); linspace(a, Comp(1, -1), Comp(28, -28));
		VecChar v(7); linspace(v, 1, 7);
		Vector<Lcomp> v1(4), v2(a.ncols());
		v1[0] = Comp(476, -476); v1[1] = Comp(504, -504);
		v1[2] = Comp(532, -532); v1[3] = Comp(560, -560);
		mul(v2, v, a);
		if (v2 != v1) SLS_ERR("failed!");
	}

	// matrix-matrix multiplication
	{
		MatComp a(7,4); linspace(a, Comp(1, -1), Comp(28, -28));
		FixCmat<Comp, 4, 7> b; her(b, a);
		CmatComp c(b.nrows(), a.ncols());
		mul(c, b, a);
		if (c(0, 0) != 3262 || c(0, 2) != 3626 || c(1, 1) != 3640 || c(1, 3) != 4032 ||
			c(2, 2) != 4046 || c(2, 3) != 4256 || c(3, 3) != 4480)
			SLS_ERR("failed!");
	}

	// v = cumsum(v)
	{
		VecInt v(4); linspace(v, 1, 4);
		VecLlong v1(v.size());
		cumsum(v1, v);
		if (v1[0] != 1 || v1[1] != 3 || v1[2] != 6 || v1[3] != 10)
			SLS_ERR("failed!");
	}
}
