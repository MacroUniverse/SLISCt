#pragma once
#include "slisc.h"

// display vectors and matrices
// don't use template so disp() can be call when debugging
// version 1
void disp(slisc::VecUchar_I &v);
void disp(slisc::VecInt_I &v);
void disp(slisc::VecDoub_I &v);
void disp(slisc::VecComp_I &v);
void disp(slisc::MatUchar_I &a);
void disp(slisc::MatInt_I &a);
void disp(slisc::MatDoub_I &a);
void disp(slisc::MatComp_I &a);
void disp(slisc::Mat3Doub_I &a);
void disp(slisc::Mat3Comp_I &a);
// version 2
void disp(slisc::VecUchar_I &v, slisc::Int_I precision);
void disp(slisc::VecInt_I &v, slisc::Int_I precision);
void disp(slisc::VecDoub_I &v, slisc::Int_I precision);
void disp(slisc::VecComp_I &v, slisc::Int_I precision);
void disp(slisc::MatUchar_I &a, slisc::Int_I precision);
void disp(slisc::MatInt_I &a, slisc::Int_I precision);
void disp(slisc::MatDoub_I &a, slisc::Int_I precision);
void disp(slisc::MatComp_I &a, slisc::Int_I precision);
void disp(slisc::Mat3Doub_I &a, slisc::Int_I precision);
void disp(slisc::Mat3Comp_I &a, slisc::Int_I precision);
// version 3
void disp(slisc::VecUchar_I &v, slisc::Long_I start, slisc::Long_I n);
void disp(slisc::VecInt_I &v, slisc::Long_I start, slisc::Long_I n);
void disp(slisc::VecDoub_I &v, slisc::Long_I start, slisc::Long_I n);
void disp(slisc::VecComp_I &v, slisc::Long_I start, slisc::Long_I n);
void disp(slisc::MatUchar_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void disp(slisc::MatInt_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void disp(slisc::MatDoub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void disp(slisc::MatComp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void disp(slisc::Mat3Doub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3);
void disp(slisc::Mat3Comp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3);
// version 4
void disp(slisc::VecUchar_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void disp(slisc::VecInt_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void disp(slisc::VecDoub_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void disp(slisc::VecComp_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void disp(slisc::MatUchar_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void disp(slisc::MatInt_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void disp(slisc::MatDoub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void disp(slisc::MatComp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void disp(slisc::Mat3Doub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3, slisc::Int_I precision);
void disp(slisc::Mat3Comp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3, slisc::Int_I precision);
