#include "../SLISC/mattsave.h"

void test_mattsave()
{
	using namespace slisc;

	// write to file
	MATTFile *pfmatt = mattOpen("test.matt", "w");

	// scalars
	Uchar s8 = 88;
	mattsave(s8, "s8", pfmatt);

	Int si = 99;
	mattsave(si, "si", pfmatt);

	Doub s(3.14159265358979323);
	mattsave(s, "s", pfmatt);

	Comp sc(s, -s);
	mattsave(sc, "sc", pfmatt);

	// vectors
	VecUchar v8(3);
	v8[0] = 1; v8[1] = 2; v8[2] = 3;
	mattsave(v8, "v8", pfmatt);
	VecInt vi(3);
	vi[0] = 1; vi[1] = 2; vi[2] = 3;
	mattsave(vi, "vi", pfmatt);
	VecDoub v(3);
	v[0] = 1.; v[1] = 2.; v[2] = 3.;
	mattsave(v, "v", pfmatt);
	VecComp vc(3);
	for (Int i = 0; i < 3; ++i)
		vc[i] = 1 + (Doub)i - I * (Doub)i;
	mattsave(vc, "vc", pfmatt);

	// matrices
	MatUchar A8(2, 3);
	A8 = 0;
	A8(0, 0) = 1; A8(0, 1) = 3; A8(0, 2) = 5; A8(1, 2) = 11;
	mattsave(A8, "A8", pfmatt);

	MatInt AI(2, 3);
	AI = 0;
	AI(0, 0) = 1; AI(0, 1) = 3; AI(0, 2) = 5; AI(1, 2) = 11;
	mattsave(AI, "AI", pfmatt);

	MatDoub A(2, 3);
	A = 0;
	A(0, 0) = 1; A(0, 1) = 3; A(0, 2) = 5; A(1, 2) = 11;
	mattsave(A, "A", pfmatt);

	MatComp C(3, 3);
	Comp *pC{ &C(0, 0) };
	for (Int i = 0; i < 9; ++i) {
		pC[i] = 1 + (Doub)i + I * (Doub)i;
	}
	mattsave(C, "C", pfmatt);

	// 3d arrays
	Mat3Doub A3;
	A3.resize(2, 2, 2);
	Doub *pA3 = &A3(0, 0, 0);
	for (Int i = 0; i < 8; ++i)
		pA3[i] = 1. + (Doub)i;
	mattsave(A3, "A3", pfmatt);

	Mat3Comp C3;
	C3.resize(2, 2, 2);
	Comp *pC3 = &C3(0, 0, 0);
	for (Int i = 0; i < 8; ++i)
		pC3[i] = Comp(1. + (Doub)i, (Doub)i);
	mattsave(C3, "C3", pfmatt);

	mattClose(pfmatt);

	// read test

	// scalars
	pfmatt = mattOpen("test.matt", "r");

	Int r_si;
	mattload(r_si, "si", pfmatt);
	if (r_si != si) error("failed!");

	Doub r_s;
	mattload(r_s, "s", pfmatt);
	if (abs(r_s-s) > 1e-15) error("failed!");

	Comp r_sc;
	mattload(r_sc, "sc", pfmatt);
	if (abs(r_sc-sc) > 1e-15) error("failed!");

	// vectors
	// TODO: Char

	VecInt r_vi;
	mattload(r_vi, "vi", pfmatt);
	if (r_vi != vi) error("failed!");

	VecDoub r_v;
	mattload(r_v, "v", pfmatt);
	r_v -= v;
	if (norm(r_v) > 1e-15) error("failed!");

	VecComp r_vc;
	mattload(r_vc, "vc", pfmatt);
	r_vc -= vc;
	if (norm(r_v) > 1e-15) error("failed!");

	// matrices
	// TODO: Char

	MatInt r_AI;
	mattload(r_AI, "AI", pfmatt);
	if (r_AI != AI) error("failed!");

	MatDoub r_A;
	mattload(r_A, "A", pfmatt);
	r_A -= A;
	if (norm(r_A) > 1e-15) error("failed!");

	MatComp r_C;
	mattload(r_C, "C", pfmatt);
	r_C -= C;
	if (norm(r_C) > 1e-15) error("failed!");

	// 3D arrays
	Mat3Doub r_A3;
	mattload(r_A3, "A3", pfmatt);
	r_A3 -= A3;
	if (norm(r_A3) > 1e-15) error("failed!");

	Mat3Comp r_C3;
	mattload(r_C3, "C3", pfmatt);
	r_C3 -= C3;
	if (norm(r_C3) > 1e-15) error("failed!");

	mattClose(pfmatt);
}
