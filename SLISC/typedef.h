#pragma once
#include <complex>
#include <string.h>

namespace slisc {

// Scalar types

typedef const int Int_I; // 32 bit integer
typedef int Int;
typedef int &Int_O, &Int_IO;
typedef const unsigned int Uint_I;
typedef unsigned int Uint;
typedef unsigned int &Uint_O, &Uint_IO;

#ifdef _MSC_VER
typedef const __int64 Llong_I; // 64 bit integer
typedef __int64 Llong;
typedef __int64 &Llong_O, &Llong_IO;
typedef const unsigned __int64 Ullong_I;
typedef unsigned __int64 Ullong;
typedef unsigned __int64 &Ullong_O, &Ullong_IO;
#else
typedef const long long int Llong_I; // 64 bit integer
typedef long long int Llong;
typedef long long int &Llong_O, &Llong_IO;
typedef const unsigned long long int Ullong_I;
typedef unsigned long long int Ullong;
typedef unsigned long long int &Ullong_O, &Ullong_IO;
#endif

#ifndef _USE_Int_AS_LONG
typedef Llong Long;
#else
typedef Int Long;
#endif
typedef const Long Long_I;
typedef Long;
typedef Long &Long_O, &Long_IO;

typedef const char Char_I; // 8 bit integer
typedef char Char;
typedef char &Char_O, &Char_IO;
typedef const unsigned char Uchar_I;
typedef unsigned char Uchar;
typedef unsigned char &Uchar_O, &Uchar_IO;

typedef const double Doub_I; // default floating type
typedef double Doub;
typedef double &Doub_O, &Doub_IO;

typedef const long double &Ldoub_I;
typedef long double Ldoub;
typedef long double &Ldoub_O, &Ldoub_IO;

typedef const std::complex<double> &Comp_I;
typedef std::complex<double> Comp;
typedef std::complex<double> &Comp_O, &Comp_IO;

typedef const bool Bool_I;
typedef bool Bool;
typedef bool &Bool_O, &Bool_IO;

// string type
typedef std::string Str;
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

}
