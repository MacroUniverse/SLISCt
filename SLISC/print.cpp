// real-time debug utilities
// these functions can be called in debugger
// DON'T MOVE INTO HEADER file or change anything unless you can call them in debuggers

#include "disp.h"
#include "slice_arith.h"
using namespace slisc;
using std::cout; using std::endl; using std::vector;

// print 1 element
Int elm(VecInt_I v, Long_I i)
{ return v(i); }

Long elm(VecLlong_I v, Long_I i)
{ return v(i); }

Doub elm(VecDoub_I v, Long_I i)
{ return v(i); }

Comp elm(VecComp_I v, Long_I i)
{ return v(i); }

Int elm(MatInt_I v, Long_I i, Long_I j)
{ return v(i, j); }

Long elm(MatLlong_I v, Long_I i, Long_I j)
{ return v(i, j); }

Doub elm(MatDoub_I v, Long_I i, Long_I j)
{ return v(i, j); }

Comp elm(MatComp_I v, Long_I i, Long_I j)
{ return v(i, j); }

Int elm(CmatInt_I v, Long_I i, Long_I j)
{ return v(i, j); }

Long elm(CmatLlong_I v, Long_I i, Long_I j)
{ return v(i, j); }

Doub elm(CmatDoub_I v, Long_I i, Long_I j)
{ return v(i, j); }

Comp elm(CmatComp_I v, Long_I i, Long_I j)
{ return v(i, j); }

Int elm(McooInt_I v, Long_I i, Long_I j)
{ return v(i, j); }

Long elm(McooLlong_I v, Long_I i, Long_I j)
{ return v(i, j); }

Doub elm(McooDoub_I v, Long_I i, Long_I j)
{ return v(i, j); }

Comp elm(McooComp_I v, Long_I i, Long_I j)
{ return v(i, j); }

Int elm(McoohInt_I v, Long_I i, Long_I j)
{ return v(i, j); }

Doub elm(McoohDoub_I v, Long_I i, Long_I j)
{ return v(i, j); }

Comp elm(McoohComp_I v, Long_I i, Long_I j)
{ return v(i, j); }

Int elm(Mat3Int_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Long elm(Mat3Llong_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Doub elm(Mat3Doub_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Comp elm(Mat3Comp_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Int elm(Cmat3Int_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Long elm(Cmat3Llong_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Doub elm(Cmat3Doub_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

Comp elm(Cmat3Comp_I v, Long_I i, Long_I j, Long_I k)
{ return v(i, j, k); }

//void elm(Mat4Int_I v) { disp(v); }
//
//void elm(Mat4Doub_I v) { disp(v); }
//
//void elm(Mat4Comp_I v) { disp(v); }

Int elm(Cmat4Int_I v, Long_I i, Long_I j, Long_I k, Long_I l)
{ return v(i, j, k, l); }

Long elm(Cmat4Llong_I v, Long_I i, Long_I j, Long_I k, Long_I l)
{ return v(i, j, k, l); }

Doub elm(Cmat4Doub_I v, Long_I i, Long_I j, Long_I k, Long_I l)
{ return v(i, j, k, l); }

Comp elm(Cmat4Comp_I v, Long_I i, Long_I j, Long_I k, Long_I l)
{ return v(i, j, k, l); }

// version 1

void print(VecInt_I v) { disp(v); }

void print(VecLlong_I v) { disp(v); }

void print(VecDoub_I v) { disp(v); }

void print(VecComp_I v) { disp(v); }

void print(MatInt_I v) { disp(v); }

void print(MatLlong_I v) { disp(v); }

void print(MatDoub_I v) { disp(v); }

void print(MatComp_I v) { disp(v); }

void print(CmatInt_I v) { disp(v); }

void print(CmatLlong_I v) { disp(v); }

void print(CmatDoub_I v) { disp(v); }

void print(CmatComp_I v) { disp(v); }

void print(Mat3Int_I v) { disp(v); }

void print(Mat3Llong_I v) { disp(v); }

void print(Mat3Doub_I v) { disp(v); }

void print(Mat3Comp_I v) { disp(v); }

void print(Cmat3Int_I v) { disp(v); }

void print(Cmat3Llong_I v) { disp(v); }

void print(Cmat3Doub_I v) { disp(v); }

void print(Cmat3Comp_I v) { disp(v); }

void print(Cmat3Imag_I v) { disp(v); }

//void print(Mat4Int_I v) { disp(v); }
//
//void print(Mat4Doub_I v) { disp(v); }
//
//void print(Mat4Comp_I v) { disp(v); }

void print(Cmat4Int_I v) { disp(v); }

void print(Cmat4Llong_I v) { disp(v); }

void print(Cmat4Doub_I v) { disp(v); }

void print(Cmat4Comp_I v) { disp(v); }



// version 2

#ifdef _MSC_VER
// VS supports calling overloaded functions when debugging
// version 2

void print(VecInt_I v, Int_I precision) { disp(v, precision); }

void print(VecLlong_I v, Int_I precision) { disp(v, precision); }

void print(VecDoub_I v, Int_I precision) { disp(v, precision); }

void print(VecComp_I v, Int_I precision) { disp(v, precision); }

void print(MatInt_I v, Int_I precision) { disp(v, precision); }

void print(MatLlong_I v, Int_I precision) { disp(v, precision); }

void print(MatDoub_I v, Int_I precision) { disp(v, precision); }

void print(MatComp_I v, Int_I precision) { disp(v, precision); }

void print(CmatInt_I v, Int_I precision) { disp(v, precision); }

void print(CmatLlong_I v, Int_I precision) { disp(v, precision); }

void print(CmatDoub_I v, Int_I precision) { disp(v, precision); }

void print(CmatComp_I v, Int_I precision) { disp(v, precision); }

void print(Mat3Int_I v, Int_I precision) { disp(v, precision); }

void print(Mat3Llong_I v, Int_I precision) { disp(v, precision); }

void print(Mat3Doub_I v, Int_I precision) { disp(v, precision); }

void print(Mat3Comp_I v, Int_I precision) { disp(v, precision); }

void print(Cmat3Int_I v, Int_I precision) { disp(v, precision); }

void print(Cmat3Llong_I v, Int_I precision) { disp(v, precision); }

void print(Cmat3Doub_I v, Int_I precision) { disp(v, precision); }

void print(Cmat3Comp_I v, Int_I precision) { disp(v, precision); }

void print(Cmat4Int_I v, Int_I precision) { disp(v, precision); }

void print(Cmat4Doub_I v, Int_I precision) { disp(v, precision); }

void print(Cmat4Comp_I v, Int_I precision) { disp(v, precision); }

// version 3

void print(VecInt_I v, Long_I start, Long_I n) { disp(slice_vec(v, start, n)); }

void print(VecLlong_I v, Long_I start, Long_I n) { disp(slice_vec(v, start, n)); }

void print(VecDoub_I v, Long_I start, Long_I n) { disp(slice_vec(v, start, n)); }

void print(VecComp_I v, Long_I start, Long_I n) { disp(slice_vec(v, start, n)); }

//void print(MatInt_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
//{ disp(slice(v, start1, start2, n1, n2)); }
//
//void print(MatDoub_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
//{ disp(v, start1, start2, n1, n2); }
//
//void print(MatComp_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
//{ disp(v, start1, start2, n1, n2); }
//
//void print(Mat3Doub_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
//{ disp(v, start1, start2, start3, n1, n2, n3); }
//
//void print(Mat3Comp_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
//{ disp(v, start1, start2, start3, n1, n2, n3); }

void print(CmatInt_I v, Long_I start1, Long_I n1, Long_I start2, Long_I n2)
{
	disp(slice(v, start1, n1, start2, n2));
}

void print(CmatLlong_I v, Long_I start1, Long_I n1, Long_I start2, Long_I n2)
{
	disp(slice(v, start1, n1, start2, n2));
}

void print(CmatDoub_I v, Long_I start1, Long_I n1, Long_I start2, Long_I n2)
{
	disp(slice(v, start1, n1, start2, n2));
}

void print(CmatComp_I v, Long_I start1, Long_I n1, Long_I start2, Long_I n2)
{
	disp(slice(v, start1, n1, start2, n2));
}

template <class Tmat3>
void print_3d(const Tmat3 &v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	if (n3 == 1) // slice12
		disp(slice(slice12(v, k), i, n1, j, n2));
	else if (n1 == 1) // slice23
		disp(slice(slice23(v, i), j, n2, k, n3));
	else // slice
		disp(slice(v, i, n1, j, n2, k, n3));
}

void print(Cmat3Int_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	print_3d(v, i, n1, j, n2, k, n3);
}

void print(Cmat3Llong_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	print_3d(v, i, n1, j, n2, k, n3);
}

void print(Cmat3Doub_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	print_3d(v, i, n1, j, n2, k, n3);
}

void print(Cmat3Comp_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	print_3d(v, i, n1, j, n2, k, n3);
}

void print(Cmat3Imag_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3)
{
	print_3d(v, i, n1, j, n2, k, n3);
}

template <class Tmat4>
void print_4d(const Tmat4 &v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3, Long_I l, Long_I n4)
{
	if (n3 == 1 && n4 == 1) // slice12
		disp(slice(slice12(v, k, l), i, n1, j, n2));
	else if (n1 == 1 && n2 == 1) // slice34
		disp(slice(slice34(v, i, j), k, n3, l, n4));
	else // slice
		disp(slice(v, i, n1, j, n2, k, n3, l, n4));
}

void print(Cmat4Int_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3, Long_I l, Long_I n4)
{
	print_4d(v, i, n1, j, n2, k, n3, l, n4);
}

void print(Cmat4Llong_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3, Long_I l, Long_I n4)
{
	print_4d(v, i, n1, j, n2, k, n3, l, n4);
}

void print(Cmat4Doub_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3, Long_I l, Long_I n4)
{
	print_4d(v, i, n1, j, n2, k, n3, l, n4);
}

void print(Cmat4Comp_I v, Long_I i, Long_I n1, Long_I j, Long_I n2, Long_I k, Long_I n3, Long_I l, Long_I n4)
{
	print_4d(v, i, n1, j, n2, k, n3, l, n4);
}

// version 4

void print(VecInt_I v, Long_I start, Long_I n, Int_I precision)
{ disp(slice_vec(v, start, n), precision); }

void print(VecLlong_I v, Long_I start, Long_I n, Int_I precision)
{ disp(slice_vec(v, start, n), precision); }

void print(VecDoub_I v, Long_I start, Long_I n, Int_I precision)
{ disp(slice_vec(v, start, n), precision); }

void print(VecComp_I v, Long_I start, Long_I n, Int_I precision)
{ disp(slice_vec(v, start, n), precision); }

//void print(MatUchar_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
//{ disp(v, start1, start2, n1, n2, precision); }
//
//void print(MatInt_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
//{ disp(v, start1, start2, n1, n2, precision); }
//
//void print(MatDoub_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
//{ disp(v, start1, start2, n1, n2, precision); }
//
//void print(MatComp_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
//{ disp(v, start1, start2, n1, n2, precision); }
//
//void print(Mat3Doub_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
//{ disp(v, start1, start2, start3, n1, n2, n3, precision); }
//
//void print(Mat3Comp_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
//{ disp(v, start1, start2, start3, n1, n2, n3, precision); }

#else // #ifdef _MSC_VER
// GCC does not fully support calling overloaded functions when debugging
// must specify number of arguments
// version 2

void print2(VecUchar_I v, Int_I precision) { disp(v, precision); }

void print2(VecInt_I v, Int_I precision) { disp(v, precision); }

void print2(VecDoub_I v, Int_I precision) { disp(v, precision); }

void print2(VecComp_I v, Int_I precision) { disp(v, precision); }

void print2(MatUchar_I v, Int_I precision) { disp(v, precision); }

void print2(MatInt_I v, Int_I precision) { disp(v, precision); }

void print2(MatDoub_I v, Int_I precision) { disp(v, precision); }

void print2(MatComp_I v, Int_I precision) { disp(v, precision); }

void print2(Mat3Doub_I v, Int_I precision) { disp(v, precision); }

void print2(Mat3Comp_I v, Int_I precision) { disp(v, precision); }

// version 3

void print3(VecUchar_I v, Long_I start, Long_I n) { disp(v, start, n); }

void print3(VecInt_I v, Long_I start, Long_I n) { disp(v, start, n); }

void print3(VecDoub_I v, Long_I start, Long_I n) { disp(v, start, n); }

void print3(VecComp_I v, Long_I start, Long_I n) { disp(v, start, n); }

void print5(MatUchar_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print5(MatInt_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print5(MatDoub_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print5(MatComp_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print7(Mat3Doub_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{ disp(v, start1, start2, start3, n1, n2, n3); }

void print7(Mat3Comp_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{ disp(v, start1, start2, start3, n1, n2, n3); }

// version 4

void print4(VecUchar_I v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print4(VecInt_I v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print4(VecDoub_I v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print4(VecComp_I v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print6(MatUchar_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print6(MatInt_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print6(MatDoub_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print6(MatComp_I v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print8(Mat3Doub_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{ disp(v, start1, start2, start3, n1, n2, n3, precision); }

void print8(Mat3Comp_I v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{ disp(v, start1, start2, start3, n1, n2, n3, precision); }

#endif // #ifdef _MSC_VER
