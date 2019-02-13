// every program that uses SLISC should include "global.h" first

#pragma once
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#define SLS_USE_MKL // use Intel MKL when possible
#ifndef NDEBUG
#define SLS_NAN_ERROR // check nan at certain places
#endif
#define SLS_FP_EXCEPT // turn on floating point exception

#ifdef SLS_USE_MKL
#include <mkl.h>
#endif
#include <limits>
#include <cmath>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace slisc {

// using std

using std::complex;
using std::vector; using std::string; using std::to_string;
using std::cout; using std::endl;
using std::ifstream; using std::ofstream;

// Scalar types

typedef int Int;
typedef const Int Int_I; // 32 bit integer
typedef Int &Int_O, &Int_IO;

typedef const unsigned int Uint_I;
typedef unsigned int Uint;
typedef unsigned int &Uint_O, &Uint_IO;

#ifdef _MSC_VER
typedef __int64 Llong;
typedef const Llong Llong_I; // 64 bit integer
typedef Llong &Llong_O, &Llong_IO;

typedef const unsigned __int64 Ullong_I;
typedef unsigned __int64 Ullong;
typedef unsigned __int64 &Ullong_O, &Ullong_IO;
#else
typedef long long int Llong;
typedef const Llong Llong_I; // 64 bit integer
typedef Llong &Llong_O, &Llong_IO;

typedef unsigned Llong Ullong;
typedef const Ullong Ullong_I;
typedef Ullong &Ullong_O, &Ullong_IO;
#endif

#ifndef SLS_USE_INT_AS_LONG
typedef Llong Long;
#else
typedef Int Long;
#endif
typedef const Long Long_I;
typedef Long &Long_O, &Long_IO;

typedef char Char;
typedef const Char Char_I; // 8 bit integer
typedef Char &Char_O, &Char_IO;

typedef unsigned char Uchar;
typedef const Uchar Uchar_I;
typedef Uchar &Uchar_O, &Uchar_IO;

typedef float Float;
typedef const Float Float_I; // default floating type
typedef Float &Float_O, &Float_IO;

typedef double Doub;
typedef const Doub Doub_I; // default floating type
typedef Doub &Doub_O, &Doub_IO;

typedef long double Ldoub;
typedef const Ldoub &Ldoub_I;
typedef Ldoub &Ldoub_O, &Ldoub_IO;

typedef complex<Float> Fcomp;
typedef const Fcomp &Fcomp_I;
typedef Fcomp &Fcomp_O, &Fcomp_IO;

typedef complex<Doub> Comp;
typedef const Comp &Comp_I;
typedef Comp &Comp_O, &Comp_IO;

typedef complex<Ldoub> Lcomp;
typedef const Lcomp &Lcomp_I;
typedef Lcomp &Lcomp_O, &Lcomp_IO;

typedef bool Bool;
typedef const Bool Bool_I;
typedef Bool &Bool_O, &Bool_IO;

// string type
typedef string Str;
typedef const Str &Str_I;
typedef Str &Str_O, &Str_IO;

// === class declarations ===
template <class T> class Vector;
template <class T> class Matrix;
template <class T> class Cmat;
template <class T, Long Nr> class FixVec;
template <class T, Long Nr, Long Nc> class FixCmat;
template <class T> class Mat3d;
template <class T> class Diag;
template <class T> class MatCoo;
template <class T> class MatCooH;
// For cuSLISC project
#ifdef _CUSLISC_
template <class T> class Gvector;
template <class T> class Gmatrix;
template <class T> class Gmat3d;
#endif

// vector and matrix types
typedef const Vector<Int> &VecInt_I;
typedef Vector<Int> VecInt;
typedef Vector<Int> &VecInt_O, &VecInt_IO;

typedef const Vector<Uint> &VecUint_I;
typedef Vector<Uint> VecUint;
typedef Vector<Uint> &VecUint_O, &VecUint_IO;

typedef const Vector<Long> &VecLong_I;
typedef Vector<Long> VecLong;
typedef Vector<Long> &VecLong_O, &VecLong_IO;

typedef const Vector<Llong> &VecLlong_I;
typedef Vector<Llong> VecLlong;
typedef Vector<Llong> &VecLlong_O, &VecLlong_IO;

typedef const Vector<Ullong> &VecUllong_I;
typedef Vector<Ullong> VecUllong;
typedef Vector<Ullong> &VecUllong_O, &VecUllong_IO;

typedef const Vector<Char> &VecChar_I;
typedef Vector<Char> VecChar;
typedef Vector<Char> &VecChar_O, &VecChar_IO;

typedef const Vector<Char*> &VecCharp_I;
typedef Vector<Char*> VecCharp;
typedef Vector<Char*> &VecCharp_O, &VecCharp_IO;

typedef const Vector<Uchar> &VecUchar_I;
typedef Vector<Uchar> VecUchar;
typedef Vector<Uchar> &VecUchar_O, &VecUchar_IO;

typedef const Vector<Doub> &VecDoub_I;
typedef Vector<Doub> VecDoub;
typedef Vector<Doub> &VecDoub_O, &VecDoub_IO;

typedef const Vector<Doub*> &VecDoubp_I;
typedef Vector<Doub*> VecDoubp;
typedef Vector<Doub*> &VecDoubp_O, &VecDoubp_IO;

typedef const Vector<Comp> &VecComp_I;
typedef Vector<Comp> VecComp;
typedef Vector<Comp> &VecComp_O, &VecComp_IO;

typedef const Vector<Bool> &VecBool_I;
typedef Vector<Bool> VecBool;
typedef Vector<Bool> &VecBool_O, &VecBool_IO;

typedef const Matrix<Int> &MatInt_I;
typedef Matrix<Int> MatInt;
typedef Matrix<Int> &MatInt_O, &MatInt_IO;

typedef const Matrix<Uint> &MatUint_I;
typedef Matrix<Uint> MatUint;
typedef Matrix<Uint> &MatUint_O, &MatUint_IO;

typedef const Matrix<Llong> &MatLlong_I;
typedef Matrix<Llong> MatLlong;
typedef Matrix<Llong> &MatLlong_O, &MatLlong_IO;

typedef const Matrix<Ullong> &MatUllong_I;
typedef Matrix<Ullong> MatUllong;
typedef Matrix<Ullong> &MatUllong_O, &MatUllong_IO;

typedef const Matrix<Char> &MatChar_I;
typedef Matrix<Char> MatChar;
typedef Matrix<Char> &MatChar_O, &MatChar_IO;

typedef const Matrix<Uchar> &MatUchar_I;
typedef Matrix<Uchar> MatUchar;
typedef Matrix<Uchar> &MatUchar_O, &MatUchar_IO;

typedef const Matrix<Doub> &MatDoub_I;
typedef Matrix<Doub> MatDoub;
typedef Matrix<Doub> &MatDoub_O, &MatDoub_IO;

typedef const Matrix<Comp> &MatComp_I;
typedef Matrix<Comp> MatComp;
typedef Matrix<Comp> &MatComp_O, &MatComp_IO;

typedef const Matrix<Bool> &MatBool_I;
typedef Matrix<Bool> MatBool;
typedef Matrix<Bool> &MatBool_O, &MatBool_IO;

typedef const Cmat<Int> &CmatInt_I;
typedef Cmat<Int> CmatInt;
typedef Cmat<Int> &CmatInt_O, &CmatInt_IO;

typedef const Cmat<Uint> &CmatUint_I;
typedef Cmat<Uint> CmatUint;
typedef Cmat<Uint> &CmatUint_O, &CmatUint_IO;

typedef const Cmat<Llong> &CmatLlong_I;
typedef Cmat<Llong> CmatLlong;
typedef Cmat<Llong> &CmatLlong_O, &CmatLlong_IO;

typedef const Cmat<Ullong> &CmatUllong_I;
typedef Cmat<Ullong> CmatUllong;
typedef Cmat<Ullong> &CmatUllong_O, &CmatUllong_IO;

typedef const Cmat<Char> &CmatChar_I;
typedef Cmat<Char> CmatChar;
typedef Cmat<Char> &CmatChar_O, &CmatChar_IO;

typedef const Cmat<Uchar> &CmatUchar_I;
typedef Cmat<Uchar> CmatUchar;
typedef Cmat<Uchar> &CmatUchar_O, &CmatUchar_IO;

typedef const Cmat<Doub> &CmatDoub_I;
typedef Cmat<Doub> CmatDoub;
typedef Cmat<Doub> &CmatDoub_O, &CmatDoub_IO;

typedef const Cmat<Comp> &CmatComp_I;
typedef Cmat<Comp> CmatComp;
typedef Cmat<Comp> &CmatComp_O, &CmatComp_IO;

typedef const Cmat<Bool> &CmatBool_I;
typedef Cmat<Bool> CmatBool;
typedef Cmat<Bool> &CmatBool_O, &CmatBool_IO;

typedef const Mat3d<Doub> &Mat3Doub_I;
typedef Mat3d<Doub> Mat3Doub;
typedef Mat3d<Doub> &Mat3Doub_O, &Mat3Doub_IO;

typedef const Mat3d<Comp> &Mat3Comp_I;
typedef Mat3d<Comp> Mat3Comp;
typedef Mat3d<Comp> &Mat3Comp_O, &Mat3Comp_IO;

// fixed-size containers

template <Long N> using FvecChar = FixVec<Char, N>;
template <Long N> using FvecChar_I = const FixVec<Char, N> &;
template <Long N> using FvecChar_O = FixVec<Char, N> &;
template <Long N> using FvecChar_IO = FixVec<Char, N> &;

template <Long N> using FvecInt = FixVec<Int, N>;
template <Long N> using FvecInt_I = const FixVec<Int, N> &;
template <Long N> using FvecInt_O = FixVec<Int, N> &;
template <Long N> using FvecInt_IO = FixVec<Int, N> &;

template <Long N> using FvecDoub = FixVec<Doub, N>;
template <Long N> using FvecDoub_I = const FixVec<Doub, N> &;
template <Long N> using FvecDoub_O = FixVec<Doub, N> &;
template <Long N> using FvecDoub_IO = FixVec<Doub, N> &;

template <Long N> using FvecComp = FixVec<Comp, N>;
template <Long N> using FvecComp_I = const FixVec<Comp, N> &;
template <Long N> using FvecComp_O = FixVec<Comp, N> &;
template <Long N> using FvecComp_IO = FixVec<Comp, N> &;

template <Long Nr, Long Nc> using FcmatChar = FixCmat<Char, Nr, Nc>;
template <Long Nr, Long Nc> using FcmatChar_I = const FixCmat<Char, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatChar_O = FixCmat<Char, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatChar_IO = FixCmat<Char, Nr, Nc> &;

template <Long Nr, Long Nc> using FcmatInt = FixCmat<Int, Nr, Nc>;
template <Long Nr, Long Nc> using FcmatInt_I = const FixCmat<Int, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatInt_O = FixCmat<Int, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatInt_IO = FixCmat<Int, Nr, Nc> &;

template <Long Nr, Long Nc> using FcmatDoub = FixCmat<Doub, Nr, Nc>;
template <Long Nr, Long Nc> using FcmatDoub_I = const FixCmat<Doub, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatDoub_O = FixCmat<Doub, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatDoub_IO = FixCmat<Doub, Nr, Nc> &;

template <Long Nr, Long Nc> using FcmatComp = FixCmat<Comp, Nr, Nc>;
template <Long Nr, Long Nc> using FcmatComp_I = const FixCmat<Comp, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatComp_O = FixCmat<Comp, Nr, Nc> &;
template <Long Nr, Long Nc> using FcmatComp_IO = FixCmat<Comp, Nr, Nc> &;

// sparse containers

typedef Diag<Int> DiagInt;
typedef const DiagInt &DiagInt_I;
typedef DiagInt &DiagInt_O, &DiagInt_IO;

typedef Diag<Doub> DiagDoub;
typedef const DiagDoub &DiagDoub_I;
typedef DiagDoub &DiagDoub_O, &DiagDoub_IO;

typedef Diag<Comp> DiagComp;
typedef const DiagComp &DiagComp_I;
typedef DiagComp &DiagComp_O, &DiagComp_IO;

typedef MatCoo<Doub> McooDoub;
typedef const McooDoub &McooDoub_I;
typedef McooDoub &McooDoub_O, &McooDoub_IO;

typedef MatCoo<Comp> McooComp;
typedef const McooComp &McooComp_I;
typedef McooComp &McooComp_O, &McooComp_IO;

typedef MatCooH<Doub> McoohDoub;
typedef const McoohDoub &McoohDoub_I;
typedef McoohDoub &McoohDoub_O, &McoohDoub_IO;

typedef MatCooH<Comp> McoohComp;
typedef const McoohComp &McoohComp_I;
typedef McoohComp &McoohComp_O, &McoohComp_IO;

// quiet NaN definition
// uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;
//Doub NaN = sqrt(-1.);
static const Doub NaN = std::numeric_limits<Doub>::quiet_NaN();

// Floating Point Exceptions for Microsoft compilers
// no exception for integer overflow
#ifdef SLS_FP_EXCEPT
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp(0, 0);
		// also: EM_INEXACT
		cw &= ~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE | EM_UNDERFLOW | EM_DENORMAL);
		_controlfp(cw, MCW_EM);
	}
};
// in case of ODR error, put this in main function;
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif
#endif

// === constants ===

const Doub PI = 3.14159265358979323;
const Doub E = 2.71828182845904524;
const Comp I(0., 1.);

// report error and pause execution
#define error(str) do{cout << "error: " << __FILE__ << ": line " << __LINE__ << ": " << str << endl; getchar();} while(0)

#define warning(str) do{cout << "warning: " << __FILE__ << ": line " << __LINE__ << ": " << str << endl;} while(0)

}
