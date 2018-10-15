// print() is same as disp(), but for calling in debugger
// DON'T MOVE INTO HEADER file or change anything unless you can call them in debuggers

#pragma once
#include "disp.h"

// version 1
void print(slisc::VecUchar_I &v);
void print(slisc::VecInt_I &v);
void print(slisc::VecDoub_I &v);
void print(slisc::VecComp_I &v);
void print(slisc::MatUchar_I &a);
void print(slisc::MatInt_I &a);
void print(slisc::MatDoub_I &a);
void print(slisc::MatComp_I &a);
void print(slisc::Mat3Doub_I &a);
void print(slisc::Mat3Comp_I &a);
// version 2
void print(slisc::VecUchar_I &v, slisc::Int_I precision);
void print(slisc::VecInt_I &v, slisc::Int_I precision);
void print(slisc::VecDoub_I &v, slisc::Int_I precision);
void print(slisc::VecComp_I &v, slisc::Int_I precision);
void print(slisc::MatUchar_I &a, slisc::Int_I precision);
void print(slisc::MatInt_I &a, slisc::Int_I precision);
void print(slisc::MatDoub_I &a, slisc::Int_I precision);
void print(slisc::MatComp_I &a, slisc::Int_I precision);
void print(slisc::Mat3Doub_I &a, slisc::Int_I precision);
void print(slisc::Mat3Comp_I &a, slisc::Int_I precision);
// version 3
void print(slisc::VecUchar_I &v, slisc::Long_I start, slisc::Long_I n);
void print(slisc::VecInt_I &v, slisc::Long_I start, slisc::Long_I n);
void print(slisc::VecDoub_I &v, slisc::Long_I start, slisc::Long_I n);
void print(slisc::VecComp_I &v, slisc::Long_I start, slisc::Long_I n);
void print(slisc::MatUchar_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void print(slisc::MatInt_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void print(slisc::MatDoub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void print(slisc::MatComp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2);
void print(slisc::Mat3Doub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3);
void print(slisc::Mat3Comp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3);
// version 4
void print(slisc::VecUchar_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void print(slisc::VecInt_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void print(slisc::VecDoub_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void print(slisc::VecComp_I &v, slisc::Long_I start, slisc::Long_I n, slisc::Int_I precision);
void print(slisc::MatUchar_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void print(slisc::MatInt_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void print(slisc::MatDoub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void print(slisc::MatComp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I n1, slisc::Long_I n2, slisc::Int_I precision);
void print(slisc::Mat3Doub_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3, slisc::Int_I precision);
void print(slisc::Mat3Comp_I &a, slisc::Long_I start1, slisc::Long_I start2, slisc::Long_I start3, slisc::Long_I n1, slisc::Long_I n2, slisc::Long_I n3, slisc::Int_I precision);
