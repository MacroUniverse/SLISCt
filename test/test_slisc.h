#pragma once
#include "../SLISC/arithmetic.h"

void test_slisc()
{
	using namespace slisc;
	// default initialize
	{
	VecDoub vDoub;
	if (vDoub.size() != 0) error("failed!");
	if (vDoub.ptr() != nullptr) error("failed!");
	MatDoub aDoub;
	if (aDoub.size() != 0) error("failed!");
	if (aDoub.nrows() != 0) error("failed!");
	if (aDoub.ncols() != 0) error("failed!");
	if (aDoub.ptr() != nullptr) error("failed!");
	Mat3Doub a3Doub;
	if (a3Doub.size() != 0) error("failed!");
	if (a3Doub.dim1() != 0) error("failed!");
	if (a3Doub.dim2() != 0) error("failed!");
	if (a3Doub.dim3() != 0) error("failed!");
	if (a3Doub.ptr() != nullptr) error("failed!");
	}

	// size initialize
	{
	VecDoub vDoub(3);
	if (vDoub.size() != 3) error("failed!");
	if (vDoub.ptr() != &vDoub[0]) error("failed!");
	MatDoub aDoub(3, 3);
	if (aDoub.size() != 9) error("failed!");
	if (aDoub.nrows() != 3) error("failed!");
	if (aDoub.ncols() != 3) error("failed!");
	if (aDoub.ptr() != &aDoub(0, 0)) error("failed!");
	Mat3Doub a3Doub(3, 3, 3);
	if (a3Doub.size() != 27) error("failed!");
	if (a3Doub.dim1() != 3) error("failed!");
	if (a3Doub.dim2() != 3) error("failed!");
	if (a3Doub.dim3() != 3) error("failed!");
	if (a3Doub.ptr() != &a3Doub(0,0,0)) error("failed!");
	}

	// const initialize
	VecDoub vDoub(3, 1.);
	if (vDoub != 1.) error("failed!");
	MatDoub aDoub(3, 3, 1.);
	if (aDoub != 1.) error("failed!");
	Mat3Doub a3Doub(3, 3, 3, 1.);
	if (a3Doub != 1.) error("failed!");

	// resize
	vDoub.resize(0);
	if (vDoub.ptr() != nullptr) error("failed!");
	vDoub.resize(4);
	if (vDoub.size() != 4) error("failed!");
	if (vDoub.ptr() != &vDoub[0]) error("failed!");
	aDoub.resize(0, 3);
	if (aDoub.ptr() != nullptr) error("failed!");
	aDoub.resize(3, 0);
	if (aDoub.ptr() != nullptr) error("failed!");
	aDoub.resize(4, 4);
	if (aDoub.size() != 16) error("failed!");
	if (aDoub.nrows() != 4) error("failed!");
	if (aDoub.ncols() != 4) error("failed!");
	if (aDoub.ptr() != &aDoub(0,0)) error("failed!");
	a3Doub.resize(0, 0, 4);
	if (a3Doub.ptr() != nullptr) error("failed!");
	a3Doub.resize(0, 4, 0);
	if (a3Doub.ptr() != nullptr) error("failed!");
	a3Doub.resize(4, 0, 0);
	if (a3Doub.ptr() != nullptr) error("failed!");
	a3Doub.resize(4, 4, 4);
	if (a3Doub.size() != 64) error("failed!");
	if (a3Doub.dim1() != 4) error("failed!");
	if (a3Doub.dim2() != 4) error("failed!");
	if (a3Doub.dim3() != 4) error("failed!");
	if (a3Doub.ptr() != &a3Doub(0,0,0)) error("failed!");

	// assignment operator
	vDoub = 1.; if (vDoub != 1.) error("failed!");
	aDoub = 1.; if (aDoub != 1.) error("failed!");
	a3Doub = 1.; if (a3Doub != 1.) error("failed!");
	VecDoub vDoub1(4);
	vDoub1 = 2.;
	vDoub = vDoub1;
	if (vDoub != vDoub1) error("failed!");
	MatDoub aDoub1(4, 4);
	aDoub1 = 2.;
	aDoub = aDoub1;
	if (aDoub != aDoub1) error("failed!");
	Mat3Doub a3Doub1(4, 4, 4);
	a3Doub1 = 2.;
	a3Doub = a3Doub1;
	if (a3Doub != a3Doub1) error("failed!");

	// move operator
	VecDoub vDoub2;
	MatDoub aDoub2;
	Mat3Doub a3Doub2;
	vDoub2 << vDoub;
	if (vDoub2 != vDoub1) error("failed!");
	if (vDoub.size() != 0) error("failed!");
	if (vDoub.ptr() != 0)  error("failed!");
	aDoub2 << aDoub;
	if (aDoub2 != aDoub1)  error("failed!");
	if (aDoub.size() != 0) error("failed!");
	if (aDoub.nrows() != 0) error("failed!");
	if (aDoub.ncols() != 0) error("failed!");
	a3Doub2 << a3Doub;
	if (a3Doub2 != a3Doub1) error("failed!");
	if (a3Doub.size() != 0) error("failed!");
	if (a3Doub.dim1() != 0) error("failed!");
	if (a3Doub.dim2() != 0) error("failed!");
	if (a3Doub.dim3() != 0) error("failed!");

	// end()
	vDoub1[vDoub1.size()-1] = 3.1;
	if (vDoub1.end() != 3.1)  error("failed!");
	if (vDoub1.end(1) != 3.1)  error("failed!");
	aDoub1(aDoub1.size()-1) = 3.1;
	if (aDoub1.end() != 3.1)  error("failed!");
	if (aDoub1.end(1) != 3.1)  error("failed!");
	a3Doub1(a3Doub1.size()-1) = 3.1;
	if (a3Doub1.end() != 3.1)  error("failed!");
	if (a3Doub1.end(1) != 3.1)  error("failed!");
	vDoub1[vDoub1.size()-2] = 3.2;
	if (vDoub1.end(2) != 3.2)  error("failed!");
	aDoub1(aDoub1.size()-2) = 3.2;
	if (aDoub1.end(2) != 3.2)  error("failed!");
	a3Doub1(a3Doub1.size()-2) = 3.2;
	if (a3Doub1.end(2) != 3.2)  error("failed!");

	// operator()
	vDoub1[vDoub1.size()-1] = 5.5;
	if ( vDoub1[vDoub1.size()-1] != 5.5 ) error("failed!");
	if (vDoub1.end() != 5.5)  error("failed!");
	aDoub1(aDoub1.nrows()-1, aDoub1.ncols()-1) = 5.5;
	if ( aDoub1(aDoub1.nrows()-1, aDoub1.ncols()-1) != 5.5 ) error("failed!");
	if (aDoub1.end() != 5.5)  error("failed!");
	a3Doub1(a3Doub1.dim1()-1, a3Doub1.dim2()-1, a3Doub1.dim3()-1) = 5.5;
	if ( a3Doub1(a3Doub1.dim1()-1, a3Doub1.dim2()-1, a3Doub1.dim3()-1) != 5.5 ) error("failed!");
	if (a3Doub1.end() != 5.5)  error("failed!");

	// operator()
	// TODO:
}
