// display a vector or matrix
// should only be called in debugger
// DON'T MOVE INTO HEADER file or change anything unless you can call them in debuggers

#include "disp.h"
using namespace slisc;
using std::cout; using std::endl; using std::vector;

// version 1

void print(VecUchar_I &v) { disp(v); }

void print(VecInt_I &v) { disp(v); }

void print(VecDoub_I &v) { disp(v); }

void print(VecComp_I &v) { disp(v); }

void print(MatUchar_I &v) { disp(v); }

void print(MatInt_I &v) { disp(v); }

void print(MatDoub_I &v) { disp(v); }

void print(MatComp_I &v) { disp(v); }

void print(Mat3Doub_I &v) { disp(v); }

void print(Mat3Comp_I &v) { disp(v); }

// version 2

void print(VecUchar_I &v, Int_I precision) { disp(v, precision); }

void print(VecInt_I &v, Int_I precision) { disp(v, precision); }

void print(VecDoub_I &v, Int_I precision) { disp(v, precision); }

void print(VecComp_I &v, Int_I precision) { disp(v, precision); }

void print(MatUchar_I &v, Int_I precision) { disp(v, precision); }

void print(MatInt_I &v, Int_I precision) { disp(v, precision); }

void print(MatDoub_I &v, Int_I precision) { disp(v, precision); }

void print(MatComp_I &v, Int_I precision) { disp(v, precision); }

void print(Mat3Doub_I &v, Int_I precision) { disp(v, precision); }

void print(Mat3Comp_I &v, Int_I precision) { disp(v, precision); }

// version 3

void print(VecUchar_I &v, Long_I start, Long_I n) { disp(v, start, n); }

void print(VecInt_I &v, Long_I start, Long_I n) { disp(v, start, n); }

void print(VecDoub_I &v, Long_I start, Long_I n) { disp(v, start, n); }

void print(VecComp_I &v, Long_I start, Long_I n) { disp(v, start, n); }

void print(MatUchar_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print(MatInt_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print(MatDoub_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print(MatComp_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{ disp(v, start1, start2, n1, n2); }

void print(Mat3Doub_I &v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{ disp(v, start1, start2, start3, n1, n2, n3); }

void print(Mat3Comp_I &v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{ disp(v, start1, start2, start3, n1, n2, n3); }

// version 4

void print(VecUchar_I &v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print(VecInt_I &v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print(VecDoub_I &v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print(VecComp_I &v, Long_I start, Long_I n, Int_I precision)
{ disp(v, start, n, precision); }

void print(MatUchar_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print(MatInt_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print(MatDoub_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print(MatComp_I &v, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{ disp(v, start1, start2, n1, n2, precision); }

void print(Mat3Doub_I &v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{ disp(v, start1, start2, start3, n1, n2, n3, precision); }

void print(Mat3Comp_I &v, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{ disp(v, start1, start2, start3, n1, n2, n3, precision); }
