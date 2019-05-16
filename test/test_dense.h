// test dense containers here
#pragma once
#include "../SLISC/arithmetic.h"

void test_dense()
{
	using namespace slisc;

	// static or constexpr members
	{
		if (VecInt::ndims() != 1) SLS_ERR("failed!");
		if (MatInt::ndims() != 2 || MatInt::major() != 'r') SLS_ERR("failed!");
		if (CmatInt::ndims() != 2 || CmatInt::major() != 'c') SLS_ERR("failed!");
		if (Mat3Doub::ndims() != 3 || Mat3Doub::major() != 'r') SLS_ERR("failed!");
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
	}

	// const initialize
	VecDoub vDoub(3, 1.);
	if (vDoub != 1.) SLS_ERR("failed!");
	MatDoub aDoub(3, 3, 1.);
	if (aDoub != 1.) SLS_ERR("failed!");
	Mat3Doub a3Doub(3, 3, 3, 1.);
	if (a3Doub != 1.) SLS_ERR("failed!");

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

	// move operator
	VecDoub vDoub2(0);
	MatDoub aDoub2(0,0);
	Mat3Doub a3Doub2(0,0,0);
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
	vDoub1[vDoub1.size()-2] = 3.2;
	if (vDoub1.end(2) != 3.2)  SLS_ERR("failed!");
	aDoub1(aDoub1.size()-2) = 3.2;
	if (aDoub1.end(2) != 3.2)  SLS_ERR("failed!");
	a3Doub1(a3Doub1.size()-2) = 3.2;
	if (a3Doub1.end(2) != 3.2)  SLS_ERR("failed!");

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

	// operator()
	// TODO:
}
