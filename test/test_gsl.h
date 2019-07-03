#pragma once
#include "../SLISC/global.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>

void test_gsl()
{
	using namespace slisc;
	Int two_ja = 1, two_ma = 1;
	Int two_jb = 1, two_mb = -1;
	Int two_jc = 0, two_mc = 0;

	// should be 1/sqrt(2) = 0.707107...
	Doub ret = gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc);
	if (abs(ret - 1. / sqrt(2)) > 1e-14)
		SLS_ERR("failed!");

	Int two_jd, two_je, two_jf, two_jg, two_jh, two_ji;
	// {ja jb jc
	//  jd je jf
	//  jg jh ji}
	two_ja = 6; two_jb = 4, two_jc = 2;
	two_jd = 8; two_je = 2; two_jf = 4;
	two_jg = 0; two_jh = 6; two_ji = 4;

	// should be 0
	ret = gsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji);
	if (abs(ret) > 1e-14)
		SLS_ERR("failed!");

	Comp in(1.5, 1.5), out;

	Int lmax = 3;
	Int n = gsl_sf_legendre_array_n(lmax);
	VecDoub legen_arr(n);
	gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, lmax, 0.6, legen_arr.ptr());
	ret = legen_arr[gsl_sf_legendre_array_index(3, 2)];
	if (abs(ret - 5.76) > 1e-14)
		SLS_ERR("failed!");

	// test legendre
	{
		Int l = 3, m = 2;
		Doub x1 = 0.4, x2 = -0.2, x3 = 0.6;
		cout << "Plm(x1) = " << gsl_sf_legendre_Plm(l, m, x1) << endl;
		cout << "Plm(x2) = " << gsl_sf_legendre_Plm(l, m, x2) << endl;
		cout << "Plm(x3) = " << gsl_sf_legendre_Plm(l, m, x3) << endl;
	}
}
