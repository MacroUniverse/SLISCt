// every program that uses SLISC should include "global.h" first

#pragma once
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define SLS_FP_EXCEPT // turn on floating point exception

#ifndef NDEBUG
#define SLS_CHECK_BOUNDS
#define SLS_CHECK_SHAPE
#endif

#include <limits>
#include <cmath>
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>

#if SLS_CPP_STD == 17
    #define SLS_CPP17
#endif

#ifdef SLS_USE_MKL
	#define MKL_Complex16 double _Complex
	#include <mkl.h>
	#define SLS_USE_CBLAS
	#define SLS_USE_LAPACKE
#else
	#ifdef SLS_USE_CBLAS
		#include <cblas.h>
	#endif
	#ifdef SLS_USE_LAPACKE
		#include <lapacke.h>
		#ifdef I // I is already defined in "/usr/include/complex.h"
			#undef I
		#endif
	#endif
#endif

namespace slisc {

// using std

using std::complex;
using std::vector; using std::string; using std::to_string;
using std::cin; using std::cout; using std::cerr; using std::endl;
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

typedef unsigned long long int Ullong;
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

typedef char32_t Char32;
typedef const Char32 Char32_I;
typedef Char32 &Char32_O, &Char32_IO;

typedef std::u32string Str32;
typedef const Str32 &Str32_I;
typedef Str32 &Str32_O, &Str32_IO;

// === class declarations ===
template <class T> class ImagNum;
template <class T> class Vector;
template <class T> class Matrix;
template <class T> class Cmat;
template <class T, Long Nr> class FixVec;
template <class T, Long Nr, Long Nc> class FixCmat;
template <class T> class Mat3d;
template <class T> class Cmat3d;
template <class T> class Cmat4d;
template <class T> class Svector;
template <class T> class Svector_c;
template <class T> class Dvector;
template <class T> class Smat;
template <class T> class Scmat;
template <class T> class Scmat_c;
template <class T> class Dmat;
template <class T> class Dcmat;
template <class T> class Jcmat;
template <class T> class Jcmat3d;
template <class T> class Jcmat4d;
template <class T> class Scmat3d;
template <class T> class Diag;
template <class T> class MatCoo;
template <class T> class MatCooH;
template <class T> class CmatObd;
template <class T> class Flm;
class Matt;

// For cuSLISC project
#ifdef _CUSLISC_
template <class T> class Gvector;
template <class T> class Gmatrix;
template <class T> class Gmat3d;
#endif

// pure imaginary number
typedef ImagNum<Float> Fimag;
typedef const Fimag Fimag_I;
typedef Fimag &Fimag_O, &Fimag_IO;

typedef ImagNum<Doub> Imag;
typedef const Imag Imag_I;
typedef Imag &Imag_O, &Imag_IO;

typedef ImagNum<Ldoub> Limag;
typedef const Limag Limag_I;
typedef Limag &Limag_O, &Limag_IO;

// vector and matrix types
typedef Vector<Int> VecInt;
typedef const VecInt &VecInt_I;
typedef VecInt &VecInt_O, &VecInt_IO;

typedef Vector<Uint> VecUint;
typedef const VecUint &VecUint_I;
typedef VecUint &VecUint_O, &VecUint_IO;

typedef Vector<Long> VecLong;
typedef const VecLong &VecLong_I;
typedef VecLong &VecLong_O, &VecLong_IO;

typedef Vector<Llong> VecLlong;
typedef const VecLlong &VecLlong_I;
typedef VecLlong &VecLlong_O, &VecLlong_IO;

typedef Vector<Ullong> VecUllong;
typedef const VecUllong &VecUllong_I;
typedef VecUllong &VecUllong_O, &VecUllong_IO;

typedef Vector<Char> VecChar;
typedef const VecChar &VecChar_I;
typedef VecChar &VecChar_O, &VecChar_IO;

typedef Vector<Char*> VecCharp;
typedef const VecCharp &VecCharp_I;
typedef VecCharp &VecCharp_O, &VecCharp_IO;

typedef Vector<Uchar> VecUchar;
typedef const VecUchar &VecUchar_I;
typedef VecUchar &VecUchar_O, &VecUchar_IO;

typedef Vector<Doub> VecDoub;
typedef const VecDoub &VecDoub_I;
typedef VecDoub &VecDoub_O, &VecDoub_IO;

typedef Vector<Doub*> VecDoubp;
typedef const VecDoubp &VecDoubp_I;
typedef VecDoubp &VecDoubp_O, &VecDoubp_IO;

typedef Vector<Imag> VecImag;
typedef const VecImag &VecImag_I;
typedef VecImag &VecImag_O, &VecImag_IO;

typedef Vector<Comp> VecComp;
typedef const VecComp &VecComp_I;
typedef VecComp &VecComp_O, &VecComp_IO;

typedef Vector<Bool> VecBool;
typedef const VecBool &VecBool_I;
typedef VecBool &VecBool_O, &VecBool_IO;

typedef Matrix<Int> MatInt;
typedef const MatInt &MatInt_I;
typedef MatInt &MatInt_O, &MatInt_IO;

typedef Matrix<Uint> MatUint;
typedef const MatUint &MatUint_I;
typedef MatUint &MatUint_O, &MatUint_IO;

typedef Matrix<Long> MatLong;
typedef const MatLong &MatLong_I;
typedef MatInt &MatLong_O, &MatLong_IO;

typedef Matrix<Llong> MatLlong;
typedef const MatLlong &MatLlong_I;
typedef MatInt &MatLlong_O, &MatLlong_IO;

typedef Matrix<Ullong> MatUllong;
typedef const MatUllong &MatUllong_I;
typedef MatUllong &MatUllong_O, &MatUllong_IO;

typedef Matrix<Char> MatChar;
typedef const MatChar &MatChar_I;
typedef MatChar &MatChar_O, &MatChar_IO;

typedef Matrix<Uchar> MatUchar;
typedef const MatUchar &MatUchar_I;
typedef MatUchar &MatUchar_O, &MatUchar_IO;

typedef Matrix<Doub> MatDoub;
typedef const MatDoub &MatDoub_I;
typedef MatDoub &MatDoub_O, &MatDoub_IO;

typedef Matrix<Imag> MatImag;
typedef const MatImag &MatImag_I;
typedef MatImag &MatImag_O, &MatImag_IO;

typedef Matrix<Comp> MatComp;
typedef const MatComp &MatComp_I;
typedef MatComp &MatComp_O, &MatComp_IO;

typedef Matrix<Bool> MatBool;
typedef const MatBool &MatBool_I;
typedef MatBool &MatBool_O, &MatBool_IO;

typedef Cmat<Int> CmatInt;
typedef const CmatInt &CmatInt_I;
typedef CmatInt &CmatInt_O, &CmatInt_IO;

typedef Cmat<Uint> CmatUint;
typedef const CmatUint &CmatUint_I;
typedef CmatUint &CmatUint_O, &CmatUint_IO;

typedef Cmat<Long> CmatLong;
typedef const CmatLong &CmatLong_I;
typedef CmatLong &CmatLong_O, &CmatLong_IO;

typedef Cmat<Llong> CmatLlong;
typedef const CmatLlong &CmatLlong_I;
typedef CmatLlong &CmatLlong_O, &CmatLlong_IO;

typedef Cmat<Ullong> CmatUllong;
typedef const CmatUllong &CmatUllong_I;
typedef CmatUllong &CmatUllong_O, &CmatUllong_IO;

typedef Cmat<Char> CmatChar;
typedef const CmatChar &CmatChar_I;
typedef CmatChar &CmatChar_O, &CmatChar_IO;

typedef Cmat<Uchar> CmatUchar;
typedef const CmatUchar &CmatUchar_I;
typedef CmatUchar &CmatUchar_O, &CmatUchar_IO;

typedef Cmat<Doub> CmatDoub;
typedef const CmatDoub &CmatDoub_I;
typedef CmatDoub &CmatDoub_O, &CmatDoub_IO;

typedef Cmat<Imag> CmatImag;
typedef const CmatImag &CmatImag_I;
typedef CmatImag &CmatImag_O, &CmatImag_IO;

typedef Cmat<Comp> CmatComp;
typedef const CmatComp &CmatComp_I;
typedef CmatComp &CmatComp_O, &CmatComp_IO;

typedef Cmat<Bool> CmatBool;
typedef const CmatBool &CmatBool_I;
typedef CmatBool &CmatBool_O, &CmatBool_IO;

typedef Mat3d<Int> Mat3Int;
typedef const Mat3Int &Mat3Int_I;
typedef Mat3Int &Mat3Int_O, &Mat3Int_IO;

typedef Mat3d<Long> Mat3Long;
typedef const Mat3Long &Mat3Long_I;
typedef Mat3Long &Mat3Long_O, &Mat3Long_IO;

typedef Mat3d<Llong> Mat3Llong;
typedef const Mat3Llong &Mat3Llong_I;
typedef Mat3Llong &Mat3Llong_O, &Mat3Llong_IO;

typedef Mat3d<Doub> Mat3Doub;
typedef const Mat3Doub &Mat3Doub_I;
typedef Mat3Doub &Mat3Doub_O, &Mat3Doub_IO;

typedef Mat3d<Comp> Mat3Comp;
typedef const Mat3Comp &Mat3Comp_I;
typedef Mat3Comp &Mat3Comp_O, &Mat3Comp_IO;

typedef Cmat3d<Int> Cmat3Int;
typedef const Cmat3Int &Cmat3Int_I;
typedef Cmat3Int &Cmat3Int_O, &Cmat3Int_IO;

typedef Cmat3d<Long> Cmat3Long;
typedef const Cmat3Long &Cmat3Long_I;
typedef Cmat3Long &Cmat3Long_O, &Cmat3Long_IO;

typedef Cmat3d<Llong> Cmat3Llong;
typedef const Cmat3Llong &Cmat3Llong_I;
typedef Cmat3Llong &Cmat3Llong_O, &Cmat3Llong_IO;

typedef Cmat3d<Llong> Cmat3Llong;
typedef const Cmat3Llong &Cmat3Llong_I;
typedef Cmat3Llong &Cmat3Llong_O, &Cmat3Llong_IO;

typedef Cmat3d<Doub> Cmat3Doub;
typedef const Cmat3Doub &Cmat3Doub_I;
typedef Cmat3Doub &Cmat3Doub_O, &Cmat3Doub_IO;

typedef Cmat3d<Comp> Cmat3Comp;
typedef const Cmat3Comp &Cmat3Comp_I;
typedef Cmat3Comp &Cmat3Comp_O, &Cmat3Comp_IO;

typedef Cmat3d<Imag> Cmat3Imag;
typedef const Cmat3Imag &Cmat3Imag_I;
typedef Cmat3Imag &Cmat3Imag_O, &Cmat3Imag_IO;

typedef Cmat4d<Int> Cmat4Int;
typedef const Cmat4Int &Cmat4Int_I;
typedef Cmat4Int &Cmat4Int_O, &Cmat4Int_IO;

typedef Cmat4d<Long> Cmat4Long;
typedef const Cmat4Long &Cmat4Long_I;
typedef Cmat4Long &Cmat4Long_O, &Cmat4Long_IO;

typedef Cmat4d<Llong> Cmat4Llong;
typedef const Cmat4Llong &Cmat4Llong_I;
typedef Cmat4Llong &Cmat4Llong_O, &Cmat4Llong_IO;

typedef Cmat4d<Doub> Cmat4Doub;
typedef const Cmat4Doub &Cmat4Doub_I;
typedef Cmat4Doub &Cmat4Doub_O, &Cmat4Doub_IO;

typedef Cmat4d<Comp> Cmat4Comp;
typedef const Cmat4Comp &Cmat4Comp_I;
typedef Cmat4Comp &Cmat4Comp_O, &Cmat4Comp_IO;

typedef Svector<Char> SvecChar;
typedef SvecChar &SvecChar_O, &SvecChar_IO;
typedef Svector_c<Char> SvecChar_c;
typedef const SvecChar_c &SvecChar_I;

typedef Svector<Int> SvecInt;
typedef const SvecInt &SvecInt_I;
typedef SvecInt &SvecInt_O, &SvecInt_IO;

typedef Svector<Long> SvecLong;
typedef const SvecLong &SvecLong_I;
typedef SvecLong &SvecLong_O, &SvecLong_IO;

typedef Svector<Llong> SvecLlong;
typedef const SvecLlong &SvecLlong_I;
typedef SvecLlong &SvecLlong_O, &SvecLlong_IO;

typedef Svector<Doub> SvecDoub;
typedef SvecDoub &SvecDoub_O, &SvecDoub_IO;
typedef Svector_c<Doub> SvecDoub_c;
typedef const SvecDoub_c &SvecDoub_I;

typedef Svector<Comp> SvecComp;
typedef SvecComp &SvecComp_O, &SvecComp_IO;
typedef Svector_c<Comp> SvecComp_c;
typedef const SvecComp_c &SvecComp_I;

typedef Dvector<Int> DvecInt;
typedef const DvecInt &DvecInt_I;
typedef DvecInt &DvecInt_O, &DvecInt_IO;

typedef Dvector<Long> DvecLong;
typedef const DvecLong &DvecLong_I;
typedef DvecLong &DvecLong_O, &DvecLong_IO;

typedef Dvector<Llong> DvecLlong;
typedef const DvecLlong &DvecLlong_I;
typedef DvecLlong &DvecLlong_O, &DvecLlong_IO;

typedef Dvector<Doub> DvecDoub;
typedef const DvecDoub &DvecDoub_I;
typedef DvecDoub &DvecDoub_O, &DvecDoub_IO;

typedef Dvector<Comp> DvecComp;
typedef const DvecComp &DvecComp_I;
typedef DvecComp &DvecComp_O, &DvecComp_IO;

typedef Smat<Int> SmatInt;
typedef const SmatInt &SmatInt_I;
typedef SmatInt &SmatInt_O, &SmatInt_IO;

typedef Smat<Long> SmatLong;
typedef const SmatLong &SmatLong_I;
typedef SmatInt &SmatLong_O, &SmatLong_IO;

typedef Smat<Llong> SmatLlong;
typedef const SmatLlong &SmatLlong_I;
typedef SmatInt &SmatLlong_O, &SmatLlong_IO;

typedef Smat<Doub> SmatDoub;
typedef const SmatDoub &SmatDoub_I;
typedef SmatDoub &SmatDoub_O, &SmatDoub_IO;

typedef Smat<Comp> SmatComp;
typedef const SmatComp &SmatComp_I;
typedef SmatComp &SmatComp_O, &SmatComp_IO;

typedef Scmat<Int> ScmatInt;
typedef ScmatInt &ScmatInt_O, &ScmatInt_IO;
typedef Scmat_c<Int> ScmatInt_c;
typedef const ScmatInt_c &ScmatInt_I;

typedef Scmat<Long> ScmatLong;
typedef ScmatInt &ScmatLong_O, &ScmatLong_IO;
typedef Scmat_c<Long> ScmatLong_c;
typedef const ScmatLong_c &ScmatLong_I;

typedef Scmat<Llong> ScmatLlong;
typedef ScmatInt &ScmatLlong_O, &ScmatLlong_IO;
typedef Scmat_c<Llong> ScmatLlong_c;
typedef const ScmatLlong_c &ScmatLlong_I;

typedef Scmat<Doub> ScmatDoub;
typedef ScmatDoub &ScmatDoub_O, &ScmatDoub_IO;
typedef Scmat_c<Doub> ScmatDoub_c;
typedef const ScmatDoub_c &ScmatDoub_I;

typedef Scmat<Comp> ScmatComp;
typedef ScmatComp &ScmatComp_O, &ScmatComp_IO;
typedef Scmat_c<Comp> ScmatComp_c;
typedef const ScmatComp_c &ScmatComp_I;

typedef Dmat<Int> DmatInt;
typedef const DmatInt &DmatInt_I;
typedef DmatInt &DmatInt_O, &DmatInt_IO;

typedef Dmat<Long> DmatLong;
typedef const DmatLong &DmatLong_I;
typedef DmatInt &DmatLong_O, &DmatLong_IO;

typedef Dmat<Llong> DmatLlong;
typedef const DmatLlong &DmatLlong_I;
typedef DmatInt &DmatLlong_O, &DmatLlong_IO;

typedef Dmat<Doub> DmatDoub;
typedef const DmatDoub &DmatDoub_I;
typedef DmatDoub &DmatDoub_O, &DmatDoub_IO;

typedef Dmat<Comp> DmatComp;
typedef const DmatComp &DmatComp_I;
typedef DmatComp &DmatComp_O, &DmatComp_IO;

typedef Dcmat<Int> DcmatInt;
typedef const DcmatInt &DcmatInt_I;
typedef DcmatInt &DcmatInt_O, &DcmatInt_IO;

typedef Dcmat<Long> DcmatLong;
typedef const DcmatLong &DcmatLong_I;
typedef DcmatInt &DcmatLong_O, &DcmatLong_IO;

typedef Dcmat<Llong> DcmatLlong;
typedef const DcmatLlong &DcmatLlong_I;
typedef DcmatInt &DcmatLlong_O, &DcmatLlong_IO;

typedef Dcmat<Doub> DcmatDoub;
typedef const DcmatDoub &DcmatDoub_I;
typedef DcmatDoub &DcmatDoub_O, &DcmatDoub_IO;

typedef Dcmat<Comp> DcmatComp;
typedef const DcmatComp &DcmatComp_I;
typedef DcmatComp &DcmatComp_O, &DcmatComp_IO;

typedef Jcmat<Int> JcmatInt;
typedef const JcmatInt &JcmatInt_I;
typedef JcmatInt &JcmatInt_O, &JcmatInt_IO;

typedef Jcmat<Long> JcmatLong;
typedef const JcmatLong &JcmatLong_I;
typedef JcmatInt &JcmatLong_O, &JcmatLong_IO;

typedef Jcmat<Llong> JcmatLlong;
typedef const JcmatLlong &JcmatLlong_I;
typedef JcmatInt &JcmatLlong_O, &JcmatLlong_IO;

typedef Jcmat<Doub> JcmatDoub;
typedef const JcmatDoub &JcmatDoub_I;
typedef JcmatDoub &JcmatDoub_O, &JcmatDoub_IO;

typedef Jcmat<Comp> JcmatComp;
typedef const JcmatComp &JcmatComp_I;
typedef JcmatComp &JcmatComp_O, &JcmatComp_IO;

typedef Jcmat3d<Int> Jcmat3Int;
typedef const Jcmat3Int &Jcmat3Int_I;
typedef Jcmat3Int &Jcmat3Int_O, &Jcmat3Int_IO;

typedef Jcmat3d<Long> Jcmat3Long;
typedef const Jcmat3Long &Jcmat3Long_I;
typedef Jcmat3Int &Jcmat3Long_O, &Jcmat3Long_IO;

typedef Jcmat3d<Llong> Jcmat3Llong;
typedef const Jcmat3Llong &Jcmat3Llong_I;
typedef Jcmat3Int &Jcmat3Llong_O, &Jcmat3Llong_IO;

typedef Jcmat3d<Doub> Jcmat3Doub;
typedef const Jcmat3Doub &Jcmat3Doub_I;
typedef Jcmat3Doub &Jcmat3Doub_O, &Jcmat3Doub_IO;

typedef Jcmat3d<Comp> Jcmat3Comp;
typedef const Jcmat3Comp &Jcmat3Comp_I;
typedef Jcmat3Comp &Jcmat3Comp_O, &Jcmat3Comp_IO;

typedef Jcmat4d<Int> Jcmat4Int;
typedef const Jcmat4Int &Jcmat4Int_I;
typedef Jcmat4Int &Jcmat4Int_O, &Jcmat4Int_IO;

typedef Jcmat4d<Long> Jcmat4Long;
typedef const Jcmat4Long &Jcmat4Long_I;
typedef Jcmat4Int &Jcmat4Long_O, &Jcmat4Long_IO;

typedef Jcmat4d<Llong> Jcmat4Llong;
typedef const Jcmat4Llong &Jcmat4Llong_I;
typedef Jcmat4Int &Jcmat4Llong_O, &Jcmat4Llong_IO;

typedef Jcmat4d<Doub> Jcmat4Doub;
typedef const Jcmat4Doub &Jcmat4Doub_I;
typedef Jcmat4Doub &Jcmat4Doub_O, &Jcmat4Doub_IO;

typedef Jcmat4d<Comp> Jcmat4Comp;
typedef const Jcmat4Comp &Jcmat4Comp_I;
typedef Jcmat4Comp &Jcmat4Comp_O, &Jcmat4Comp_IO;

typedef Scmat3d<Int> Scmat3Int;
typedef const Scmat3Int &Scmat3Int_I;
typedef Scmat3Int &Scmat3Int_O, &Scmat3Int_IO;

typedef Scmat3d<Long> Scmat3Long;
typedef const Scmat3Long &Scmat3Long_I;
typedef Scmat3Long &Scmat3Long_O, &Scmat3Long_IO;

typedef Scmat3d<Llong> Scmat3Llong;
typedef const Scmat3Llong &Scmat3Llong_I;
typedef Scmat3Llong &Scmat3Llong_O, &Scmat3Llong_IO;

typedef Scmat3d<Doub> Scmat3Doub;
typedef const Scmat3Doub &Scmat3Doub_I;
typedef Scmat3Doub &Scmat3Doub_O, &Scmat3Doub_IO;

typedef Scmat3d<Comp> Scmat3Comp;
typedef const Scmat3Comp &Scmat3Comp_I;
typedef Scmat3Comp &Scmat3Comp_O, &Scmat3Comp_IO;

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

typedef MatCoo<Int> McooInt;
typedef const McooInt &McooInt_I;
typedef McooInt &McooInt_O, &McooInt_IO;

typedef MatCoo<Long> McooLong;
typedef const McooLong &McooLong_I;
typedef McooLong &McooLong_O, &McooLong_IO;

typedef MatCoo<Llong> McooLlong;
typedef const McooLlong &McooLlong_I;
typedef McooLlong &McooLlong_O, &McooLlong_IO;

typedef MatCoo<Doub> McooDoub;
typedef const McooDoub &McooDoub_I;
typedef McooDoub &McooDoub_O, &McooDoub_IO;

typedef MatCoo<Imag> McooImag;
typedef const McooImag &McooImag_I;
typedef McooImag &McooImag_O, &McooImag_IO;

typedef MatCoo<Comp> McooComp;
typedef const McooComp &McooComp_I;
typedef McooComp &McooComp_O, &McooComp_IO;

typedef MatCooH<Int> McoohInt;
typedef const McoohInt &McoohInt_I;
typedef McoohInt &McoohInt_O, &McoohInt_IO;

typedef MatCooH<Doub> McoohDoub;
typedef const McoohDoub &McoohDoub_I;
typedef McoohDoub &McoohDoub_O, &McoohDoub_IO;

typedef MatCooH<Comp> McoohComp;
typedef const McoohComp &McoohComp_I;
typedef McoohComp &McoohComp_O, &McoohComp_IO;

typedef CmatObd<Int> CmobdInt;
typedef const CmobdInt &CmobdInt_I;
typedef CmobdInt &CmobdInt_O, &CmobdInt_IO;

typedef CmatObd<Doub> CmobdDoub;
typedef const CmobdDoub &CmobdDoub_I;
typedef CmobdDoub &CmobdDoub_O, &CmobdDoub_IO;

typedef CmatObd<Imag> CmobdImag;
typedef const CmobdImag &CmobdImag_I;
typedef CmobdImag &CmobdImag_O, &CmobdImag_IO;

typedef CmatObd<Comp> CmobdComp;
typedef const CmobdComp &CmobdComp_I;
typedef CmobdComp &CmobdComp_O, &CmobdComp_IO;

typedef Flm<Doub> FlmDoub;
typedef const FlmDoub &FlmDoub_I;
typedef FlmDoub &FlmDoub_O, &FlmDoub_IO;

typedef Flm<Comp> FlmComp;
typedef const FlmComp &FlmComp_I;
typedef FlmComp &FlmComp_O, &FlmComp_IO;

typedef const Matt &Matt_I;
typedef Matt &Matt_O, &Matt_IO;

template <class T> using vector_I = const vector<T> &;
template <class T> using vector_O = vector<T> &;
template <class T> using vector_IO = vector<T> &;

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
        // also: EM_INEXACT, EM_UNDERFLOW
        cw &= ~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL);
        _controlfp(cw, MCW_EM);
    }
};
// in case of ODR error, put this in main function;
// turn_on_floating_exceptions yes_turn_on_floating_exceptions;turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif
#endif

// === constants ===

constexpr Doub PI = 3.14159265358979323;
constexpr Doub E = 2.71828182845904524;

// report error and pause execution
void pause(Doub_I t);
#ifndef SLS_ERR
#define SLS_ERR(str) do{cout << "error: " << __FILE__ << ": line " << __LINE__ << ": " << str << endl; pause(20);} while(0)
#endif

#define SLS_WARN(str) do{cout << "warning: " << __FILE__ << ": line " << __LINE__ << ": " << str << endl;} while(0)

} // namespace slisc
