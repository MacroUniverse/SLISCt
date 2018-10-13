#pragma once
#include "slisc.h"

// display vectors and matrices
// don't use template so disp() can be call when debugging
// version 1
void disp(VecUchar_I &v);
void disp(VecInt_I &v);
void disp(VecDoub_I &v);
void disp(VecComp_I &v);
void disp(MatUchar_I &a);
void disp(MatInt_I &a);
void disp(MatDoub_I &a);
void disp(MatComp_I &a);
void disp(Mat3Doub_I &a);
void disp(Mat3Comp_I &a);
// version 2
void disp(VecUchar_I &v, Int_I precision);
void disp(VecInt_I &v, Int_I precision);
void disp(VecDoub_I &v, Int_I precision);
void disp(VecComp_I &v, Int_I precision);
void disp(MatUchar_I &a, Int_I precision);
void disp(MatInt_I &a, Int_I precision);
void disp(MatDoub_I &a, Int_I precision);
void disp(MatComp_I &a, Int_I precision);
void disp(Mat3Doub_I &a, Int_I precision);
void disp(Mat3Comp_I &a, Int_I precision);
// version 3
void disp(VecUchar_I &v, Long_I start, Long_I n);
void disp(VecInt_I &v, Long_I start, Long_I n);
void disp(VecDoub_I &v, Long_I start, Long_I n);
void disp(VecComp_I &v, Long_I start, Long_I n);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
// version 4
void disp(VecUchar_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecInt_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecDoub_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecComp_I &v, Long_I start, Long_I n, Int_I precision);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);