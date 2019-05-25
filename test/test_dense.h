// test dense containers here
#pragma once
#include "../SLISC/arithmetic.h"

void test_dense()
{
	using namespace slisc;

	// static or constexpr members
	{
		if (ndims<VecInt>() != 1) SLS_ERR("failed!");
		if (ndims<MatInt>() != 2 || major<MatInt>() != 'r') SLS_ERR("failed!");
		if (ndims<CmatInt>() != 2 || major<CmatInt>() != 'c') SLS_ERR("failed!");
		if (ndims<Mat3Doub>() != 3 || major<Mat3Doub>() != 'r') SLS_ERR("failed!");
		if (ndims<Cmat3Doub>() != 3 || major<Cmat3Doub>() != 'c') SLS_ERR("failed!");
	}

	// size initialize
	{
	VecDoub vDoub(3);
	if (vDoub.size() != 3) SLS_ERR("failed!");
	if (vDoub.ptr() != &vDoub[0]) SLS_ERR("failed!");
	MatDoub aDoub(3, 3);
	if (aDoub.size() != 9) SLS_ERR("failed!");
	if (aDoub.nrows() != 3) SLS_ERR("failed!");
	if (aDoub.ncols() != 3) SLS_ERR("failed!");
	if (aDoub.ptr() != &aDoub(0, 0)) SLS_ERR("failed!");
	Mat3Doub a3Doub(3, 3, 3);
	if (a3Doub.size() != 27) SLS_ERR("failed!");
	if (a3Doub.dim1() != 3) SLS_ERR("failed!");
	if (a3Doub.dim2() != 3) SLS_ERR("failed!");
	if (a3Doub.dim3() != 3) SLS_ERR("failed!");
	if (a3Doub.ptr() != &a3Doub(0,0,0)) SLS_ERR("failed!");
	Cmat3Doub c3Doub(3, 3, 3);
	if (c3Doub.size() != 27) SLS_ERR("failed!");
	if (c3Doub.dim1() != 3) SLS_ERR("failed!");
	if (c3Doub.dim2() != 3) SLS_ERR("failed!");
	if (c3Doub.dim3() != 3) SLS_ERR("failed!");
	if (c3Doub.ptr() != &c3Doub(0, 0, 0)) SLS_ERR("failed!");
	}

	// const initialize
	VecDoub vDoub(3, 1.);
	if (vDoub != 1.) SLS_ERR("failed!");
	MatDoub aDoub(3, 3, 1.);
	if (aDoub != 1.) SLS_ERR("failed!");
	Mat3Doub a3Doub(3, 3, 3, 1.);
	if (a3Doub != 1.) SLS_ERR("failed!");
	Cmat3Doub c3Doub(3, 3, 3, 1.);
	if (c3Doub != 1.) SLS_ERR("failed!");

	// resize
	vDoub.resize(0);
	if (vDoub.size() != 0) SLS_ERR("failed!");
	vDoub.resize(4);
	if (vDoub.size() != 4) SLS_ERR("failed!");
	if (vDoub.ptr() != &vDoub[0]) SLS_ERR("failed!");
	aDoub.resize(0, 3);
	if (aDoub.size() != 0) SLS_ERR("failed!");
	aDoub.resize(3, 0);
	if (aDoub.size() != 0) SLS_ERR("failed!");
	aDoub.resize(4, 4);
	if (aDoub.size() != 16) SLS_ERR("failed!");
	if (aDoub.nrows() != 4) SLS_ERR("failed!");
	if (aDoub.ncols() != 4) SLS_ERR("failed!");
	if (aDoub.ptr() != &aDoub(0,0)) SLS_ERR("failed!");
	a3Doub.resize(0, 0, 4);
	if (a3Doub.size() != 0) SLS_ERR("failed!");
	a3Doub.resize(0, 4, 0);
	if (a3Doub.size() != 0) SLS_ERR("failed!");
	a3Doub.resize(4, 0, 0);
	if (a3Doub.size() != 0) SLS_ERR("failed!");
	a3Doub.resize(4, 4, 4);
	if (a3Doub.size() != 64) SLS_ERR("failed!");
	if (a3Doub.dim1() != 4) SLS_ERR("failed!");
	if (a3Doub.dim2() != 4) SLS_ERR("failed!");
	if (a3Doub.dim3() != 4) SLS_ERR("failed!");
	if (a3Doub.ptr() != &a3Doub(0,0,0)) SLS_ERR("failed!");
	c3Doub.resize(0, 0, 4);
	if (c3Doub.size() != 0) SLS_ERR("failed!");
	c3Doub.resize(0, 4, 0);
	if (c3Doub.size() != 0) SLS_ERR("failed!");
	c3Doub.resize(4, 0, 0);
	if (c3Doub.size() != 0) SLS_ERR("failed!");
	c3Doub.resize(4, 4, 4);
	if (c3Doub.size() != 64) SLS_ERR("failed!");
	if (c3Doub.dim1() != 4) SLS_ERR("failed!");
	if (c3Doub.dim2() != 4) SLS_ERR("failed!");
	if (c3Doub.dim3() != 4) SLS_ERR("failed!");
	if (c3Doub.ptr() != &c3Doub(0, 0, 0)) SLS_ERR("failed!");

	// resize and copy old data
	{
		Vbase<Int> v(3);
		v(0) = 100; v(1) = 101; v(2) = 102;
		v.resize_cpy(4);
		if (v(0) != 100 || v(1) != 101 || v(2) != 102)
			SLS_ERR("failed!");
		v.resize_cpy(2);
		if (v(0) != 100 || v(1) != 101)
			SLS_ERR("failed!");
	}

	// assignment operator
	vDoub = 1.; if (vDoub != 1.) SLS_ERR("failed!");
	aDoub = 1.; if (aDoub != 1.) SLS_ERR("failed!");
	a3Doub = 1.; if (a3Doub != 1.) SLS_ERR("failed!");
	VecDoub vDoub1(4);
	vDoub1 = 2.;
	vDoub = vDoub1;
	if (vDoub != vDoub1) SLS_ERR("failed!");
	MatDoub aDoub1(4, 4);
	aDoub1 = 2.;
	aDoub = aDoub1;
	if (aDoub != aDoub1) SLS_ERR("failed!");
	Mat3Doub a3Doub1(4, 4, 4);
	a3Doub1 = 2.;
	a3Doub = a3Doub1;
	if (a3Doub != a3Doub1) SLS_ERR("failed!");
	Cmat3Doub c3Doub1(4, 4, 4);
	c3Doub1 = 2.;
	c3Doub = c3Doub1;
	if (c3Doub != c3Doub1) SLS_ERR("failed!");

	// move operator
	VecDoub vDoub2(0);
	MatDoub aDoub2(0,0);
	Mat3Doub a3Doub2(0,0,0);
	Cmat3Doub c3Doub2(0, 0, 0);
	vDoub2 << vDoub;
	if (vDoub2 != vDoub1) SLS_ERR("failed!");
	if (vDoub.size() != 0) SLS_ERR("failed!");
	aDoub2 << aDoub;
	if (aDoub2 != aDoub1)  SLS_ERR("failed!");
	if (aDoub.size() != 0) SLS_ERR("failed!");
	if (aDoub.nrows() != 0) SLS_ERR("failed!");
	if (aDoub.ncols() != 0) SLS_ERR("failed!");
	a3Doub2 << a3Doub;
	if (a3Doub2 != a3Doub1) SLS_ERR("failed!");
	if (a3Doub.size() != 0) SLS_ERR("failed!");
	if (a3Doub.dim1() != 0) SLS_ERR("failed!");
	if (a3Doub.dim2() != 0) SLS_ERR("failed!");
	if (a3Doub.dim3() != 0) SLS_ERR("failed!");
	c3Doub2 << c3Doub;
	if (c3Doub2 != c3Doub1) SLS_ERR("failed!");
	if (a3Doub.size() != 0) SLS_ERR("failed!");
	if (a3Doub.dim1() != 0) SLS_ERR("failed!");
	if (a3Doub.dim2() != 0) SLS_ERR("failed!");
	if (a3Doub.dim3() != 0) SLS_ERR("failed!");

	// end()
	vDoub1[vDoub1.size()-1] = 3.1;
	if (vDoub1.end() != 3.1)  SLS_ERR("failed!");
	if (vDoub1.end(1) != 3.1)  SLS_ERR("failed!");
	aDoub1(aDoub1.size()-1) = 3.1;
	if (aDoub1.end() != 3.1)  SLS_ERR("failed!");
	if (aDoub1.end(1) != 3.1)  SLS_ERR("failed!");
	a3Doub1(a3Doub1.size()-1) = 3.1;
	if (a3Doub1.end() != 3.1)  SLS_ERR("failed!");
	if (a3Doub1.end(1) != 3.1)  SLS_ERR("failed!");
	c3Doub1(c3Doub1.size() - 1) = 3.1;
	if (c3Doub1.end() != 3.1)  SLS_ERR("failed!");
	if (c3Doub1.end(1) != 3.1)  SLS_ERR("failed!");
	vDoub1[vDoub1.size()-2] = 3.2;
	if (vDoub1.end(2) != 3.2)  SLS_ERR("failed!");
	aDoub1(aDoub1.size()-2) = 3.2;
	if (aDoub1.end(2) != 3.2)  SLS_ERR("failed!");
	a3Doub1(a3Doub1.size()-2) = 3.2;
	if (a3Doub1.end(2) != 3.2)  SLS_ERR("failed!");
	c3Doub1(c3Doub1.size() - 2) = 3.2;
	if (c3Doub1.end(2) != 3.2)  SLS_ERR("failed!");

	// operator()
	vDoub1[vDoub1.size()-1] = 5.5;
	if ( vDoub1[vDoub1.size()-1] != 5.5 ) SLS_ERR("failed!");
	if (vDoub1.end() != 5.5)  SLS_ERR("failed!");
	aDoub1(aDoub1.nrows()-1, aDoub1.ncols()-1) = 5.5;
	if ( aDoub1(aDoub1.nrows()-1, aDoub1.ncols()-1) != 5.5 ) SLS_ERR("failed!");
	if (aDoub1.end() != 5.5)  SLS_ERR("failed!");
	a3Doub1(a3Doub1.dim1()-1, a3Doub1.dim2()-1, a3Doub1.dim3()-1) = 5.5;
	if ( a3Doub1(a3Doub1.dim1()-1, a3Doub1.dim2()-1, a3Doub1.dim3()-1) != 5.5 ) SLS_ERR("failed!");
	if (a3Doub1.end() != 5.5)  SLS_ERR("failed!");

	Cmat3Doub c3d(2, 3, 4);
	linspace(c3d, 1, 24);
	if (c3d(0, 0, 0) != 1 || c3d(1,0,0) != 2)
		SLS_ERR("failed!");
	if (c3d(0, 1, 0) != 3 || c3d(0, 2, 0) != 5)
		SLS_ERR("failed!");
	if (c3d(0, 0, 1) != 7 || c3d(0, 0, 2) != 13)
		SLS_ERR("failed!");

	// operator()
	// TODO:
}
