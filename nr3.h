// type and class definitions and dependencies
// this header file can be used alone

#pragma once

#ifdef _DEBUG
// this will not check the last index
#define _CHECKBOUNDS_ 1
#endif

// all the system #include's we'll ever need
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

// Scalar types

typedef const int Int_I; // 32 bit integer
typedef int Int, Int_O, Int_IO;
typedef const unsigned int Uint_I;
typedef unsigned int Uint, Uint_O, Uint_IO;

#ifdef _MSC_VER
typedef const __int64 Llong_I; // 64 bit integer
typedef __int64 Llong, Llong_O, Llong_IO;
typedef const unsigned __int64 Ullong_I;
typedef unsigned __int64 Ullong, Ullong_O, Ullong_IO;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

#ifndef _USE_Int_AS_LONG
typedef Llong Long;
#else
typedef Int Long;
#endif
typedef const Long Long_I;
typedef Long Long_O, Long_IO;

typedef const char Char_I; // 8 bit integer
typedef char Char, Char_O, CharIO;
typedef const unsigned char Uchar_I;
typedef unsigned char Uchar, Uchar_O, Uchar_IO;

typedef const double Doub_I; // default floating type
typedef double Doub, Doub_O, Doub_IO;
typedef const long double Ldoub_I;
typedef long double Ldoub, Ldoub_O, Ldoub_IO;

typedef const std::complex<double> Comp_I;
typedef std::complex<double> Comp, Comp_O, Comp_IO;

typedef const bool Bool_I;
typedef bool Bool, Bool_O, Bool_IO;

// NaN definition
static const Doub NaN = std::numeric_limits<Doub>::quiet_NaN();

// report error and pause execution
#define error(str) {std::cout << "error: " << __FILE__ << ": line " << __LINE__ << ": " << str << std::endl; getchar();}

// macro-like inline functions
template<class T>
inline void nrmemset(T *dest, const T val, Long_I n)
{
	T *end = dest + n;
	for (; dest < end; ++dest)
		*dest = val;
}

// Vector Class

template <class T>
class NRvector {
private:
	Long nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(Long_I n);
	NRvector(Long_I n, const T &a);	//initialize to constant value
	NRvector(Long_I n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor forbidden
	inline NRvector & operator=(const NRvector &rhs);	// copy assignment
	inline NRvector & operator=(const T rhs);  // assign to constant value
	inline void operator<<(NRvector &rhs); // move data and rhs.resize(0)
	inline T & operator[](Long_I i);	//i'th element
	inline const T & operator[](Long_I i) const;
	inline Long_I size() const;
	inline void resize(Long_I newn); // resize (contents not preserved)
	template <class T1>
	inline void resize(const NRvector<T1> &v);
	~NRvector();
};

template <class T>
NRvector<T>::NRvector() : nn(0), v(nullptr) {}

template <class T>
NRvector<T>::NRvector(Long_I n) : nn(n), v(n>0 ? new T[n] : nullptr) {}

template <class T>
NRvector<T>::NRvector(Long_I n, const T& a) : NRvector(n)
{
	nrmemset(v, a, n);
}

template <class T>
NRvector<T>::NRvector(Long_I n, const T *a) : NRvector(n)
{
	memcpy(v, a, n*sizeof(T));
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!")
}

template <class T>
inline NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!")
	resize(rhs.nn);
	memcpy(v, rhs.v, nn*sizeof(T));
	return *this;
}

template <class T>
inline NRvector<T> & NRvector<T>::operator=(const T rhs)
{
	if (nn) nrmemset(v, rhs, nn);
	return *this;
}

template <class T>
inline void NRvector<T>::operator<<(NRvector<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!")
	if (v != nullptr) delete[] v;
	nn = rhs.nn;
	v = rhs.v;
	rhs.v = nullptr;
	rhs.nn = 0;
}

template <class T>
inline T & NRvector<T>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn)
	error("NRvector subscript out of bounds")
#endif
	return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn)
	error("NRvector subscript out of bounds")
#endif
	return v[i];
}

template <class T>
inline Long_I NRvector<T>::size() const
{
	return nn;
}

template <class T>
inline void NRvector<T>::resize(Long_I newn)
{
	if (newn != nn) {
		if (v != nullptr) delete[] v;
		nn = newn;
		v = nn > 0 ? new T[nn] : nullptr;
	}
}

template<class T>
template<class T1>
inline void NRvector<T>::resize(const NRvector<T1>& v)
{
	resize(v.size());
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != nullptr) delete[] v;
}

// Matrix Class

template <class T>
class NRmatrix {
private:
	Long nn;
	Long mm;
	T **v;
	inline T **data_alloc(Long_I n, Long_I m);
	inline void data_free();
public:
	NRmatrix();
	NRmatrix(Long_I n, Long_I m);
	NRmatrix(Long_I n, Long_I m, const T &a);	//Initialize to constant
	NRmatrix(Long_I n, Long_I m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	inline NRmatrix & operator=(const NRmatrix &rhs);	//assignment
	inline NRmatrix & operator=(const T rhs);
	inline void operator<<(NRmatrix &rhs); // move data and rhs.resize(0, 0)
	inline T* operator[](Long_I i);	//subscripting: pointer to row i
	inline const T* operator[](Long_I i) const;
	inline Long nrows() const;
	inline Long ncols() const;
	inline void resize(Long_I newn, Long_I newm); // resize (contents not preserved)
	template <class T1>
	inline void resize(const NRmatrix<T1> &a);
	~NRmatrix();
};

template <class T>
inline T ** NRmatrix<T>::data_alloc(Long_I n, Long_I m)
{
	if (m*n == 0) return nullptr;
	T **p = new T*[n];
	p[0] = new T[m*n];
	for (Long i = 1; i<n; i++)
		p[i] = p[i-1] + m;
	return p;
}

template <class T>
inline void NRmatrix<T>::data_free()
{
	if (v != nullptr) {
		delete v[0]; delete v;
	}
}

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(nullptr) {}

template <class T>
NRmatrix<T>::NRmatrix(Long_I n, Long_I m) : nn(n), mm(m), v(data_alloc(n, m)) {}

template <class T>
NRmatrix<T>::NRmatrix(Long_I n, Long_I m, const T &a) : NRmatrix(n, m)
{
	nrmemset(v[0], a, n*m);
}

template <class T>
NRmatrix<T>::NRmatrix(Long_I n, Long_I m, const T *a) : NRmatrix(n, m)
{
	memcpy(v[0], a, n*m*sizeof(T));
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!")
}

template <class T>
inline NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!")
	resize(rhs.nn, rhs.mm);
	memcpy(v[0], rhs.v[0], nn*mm*sizeof(T));
	return *this;
}

template <class T>
inline NRmatrix<T> & NRmatrix<T>::operator=(const T rhs)
{
	Long N = nn*mm;
	if (N) nrmemset(v[0], rhs, N);
	return *this;
}

template <class T>
inline void NRmatrix<T>::operator<<(NRmatrix<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!")
	data_free();
	nn = rhs.nn; mm = rhs.mm; v = rhs.v;
	rhs.v = nullptr;
	rhs.nn = rhs.mm = 0;
}

template <class T>
inline T* NRmatrix<T>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn)
	error("NRmatrix subscript out of bounds")
#endif
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn)
	error("NRmatrix subscript out of bounds")
#endif
	return v[i];
}

template <class T>
inline Long NRmatrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline Long NRmatrix<T>::ncols() const
{
	return mm;
}

template <class T>
inline void NRmatrix<T>::resize(Long_I newn, Long_I newm)
{
	if (newn != nn || newm != mm) {
		data_free();
		nn = newn;
		mm = newm;
		v = data_alloc(nn, mm);
	}
}

template <class T>
template <class T1>
inline void NRmatrix<T>::resize(const NRmatrix<T1> &a)
{
	resize(a.nrows(), a.ncols());
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	data_free();
}

// 3D Matrix Class

template <class T>
class NRmat3d {
private:
	Long nn;
	Long mm;
	Long kk;
	T ***v;
	inline T *** data_alloc(Long_I n, Long_I m, Long_I k);
	inline void data_free();
public:
	NRmat3d();
	NRmat3d(Long_I n, Long_I m, Long_I k);
	NRmat3d(Long_I n, Long_I m, Long_I k, const T a);
	NRmat3d(const NRmat3d &rhs);   // Copy constructor
	inline NRmat3d & operator=(const NRmat3d &rhs);	//assignment
	inline NRmat3d & operator=(const T rhs);
	inline void operator<<(NRmat3d &rhs); // move data and rhs.resize(0, 0, 0)
	inline void resize(Long_I n, Long_I m, Long_I k);
	template <class T1>
	inline void resize(const NRmat3d<T1> &a);
	inline T** operator[](Long_I i);	//subscripting: pointer to row i
	inline const T* const * operator[](Long_I i) const;
	inline Long dim1() const;
	inline Long dim2() const;
	inline Long dim3() const;
	~NRmat3d();
};

template <class T>
inline T *** NRmat3d<T>::data_alloc(Long_I n, Long_I m, Long_I k)
{
	if (n*m*k == 0) return nullptr;
	Long i;
	Long mn = m*n;
	T ***p = new T**[n];
	T **p0 = p[0] = new T*[mn];
	p[0][0] = new T[mn*k];
	for(i = 1; i < n; ++i)
		p[i] = p[i-1] + m;
	for(i = 1; i < mn; ++i)
		p0[i] = p0[i-1] + k;
	return p;
}

template <class T>
inline void NRmat3d<T>::data_free()
{
	if (v != nullptr) {
		delete v[0][0]; delete v[0]; delete v;
	}
}

template <class T>
NRmat3d<T>::NRmat3d(): nn(0), mm(0), kk(0), v(nullptr) {}

template <class T>
NRmat3d<T>::NRmat3d(Long_I n, Long_I m, Long_I k) : nn(n), mm(m), kk(k), v(data_alloc(n, m, k)) {}

template <class T>
NRmat3d<T>::NRmat3d(Long_I n, Long_I m, Long_I k, const T a) : NRmat3d(n, m, k)
{
	nrmemset(v, a, n*m*k);
}

template <class T>
NRmat3d<T>::NRmat3d(const NRmat3d<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!")
}

template <class T>
inline NRmat3d<T> &NRmat3d<T>::operator=(const NRmat3d<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!")
	resize(rhs.nn, rhs.mm, rhs.kk);
	memcpy(v[0][0], rhs.v[0][0], nn*mm*kk*sizeof(T));
	return *this;
}

template <class T>
inline NRmat3d<T> & NRmat3d<T>::operator=(const T rhs)
{
	Long N = nn*mm*kk;
	if (N) nrmemset(v[0][0], rhs, N);
	return *this;
}

template <class T>
inline void NRmat3d<T>::operator<<(NRmat3d<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!")
	data_free();
	nn = rhs.nn; mm = rhs.mm; kk = rhs.kk; v = rhs.v;
	rhs.v = nullptr;
	rhs.nn = rhs.mm = rhs.kk = 0;
}

template <class T>
inline void NRmat3d<T>::resize(Long_I n, Long_I m, Long_I k)
{
	if (n != nn || m != mm || k != kk) {
		data_free();
		nn = n; mm = m; kk = k;
		v = data_alloc(n, m, k);
	}
}

template <class T>
template <class T1>
inline void NRmat3d<T>::resize(const NRmat3d<T1> &a) { resize(a.dim1(), a.dim2(), a.dim3()); }

template <class T>
inline T** NRmat3d<T>::operator[](Long_I i) { return v[i]; }

template <class T>
inline const T* const * NRmat3d<T>::operator[](Long_I i) const { return v[i]; }

template <class T>
inline Long NRmat3d<T>::dim1() const { return nn; }

template <class T>
inline Long NRmat3d<T>::dim2() const { return mm; }

template <class T>
inline Long NRmat3d<T>::dim3() const { return kk; }

template <class T>
NRmat3d<T>::~NRmat3d() { data_free(); }

// Matric and vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Comp> VecComp_I;
typedef NRvector<Comp> VecComp, VecComp_O, VecComp_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Comp> MatComp_I;
typedef NRmatrix<Comp> MatComp, MatComp_O, MatComp_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

typedef const NRmat3d<Doub> Mat3Doub_I;
typedef NRmat3d<Doub> Mat3Doub, Mat3Doub_O, Mat3Doub_IO;

typedef const NRmat3d<Comp> Mat3Comp_I;
typedef NRmat3d<Comp> Mat3Comp, Mat3Comp_O, Mat3Comp_IO;

// macro-like functions (don't use them in your code ever, write similar utilities in "nr3plus.h")

template<class T>
inline T SQR(const T a) { return a*a; }

template<class T>
inline const T &MAX(const T &a, const T &b)
{ return b > a ? (b) : (a); }

inline float MAX(const double &a, const float &b)
{ return b > a ? (b) : float(a); }

inline float MAX(const float &a, const double &b)
{ return b > a ? float(b) : (a); }

template<class T>
inline const T &MIN(const T &a, const T &b)
{ return b < a ? (b) : (a); }

inline float MIN(const double &a, const float &b)
{ return b < a ? (b) : float(a); }

inline float MIN(const float &a, const double &b)
{ return b < a ? float(b) : (a); }

template<class T>
inline T SIGN(const T &a, const T &b)
{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

inline float SIGN(const float &a, const double &b)
{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

inline float SIGN(const double &a, const float &b)
{ return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a)); }

template<class T>
inline void SWAP(T &a, T &b)
{ T dum = a; a = b; b = dum; }

// NR3 compatibility issues

struct Complex
{ Complex() { error("Replace all \"Complex\" with \"Comp!\"") } };

struct VecComplex
{ VecComplex() { error("Replace all \"Complex\" with \"Comp!\"") } };

struct MatComplex
{ MatComplex() { error("Replace all \"Complex\" with \"Comp!\"") } };

struct Mat3DDoub
{ Mat3DDoub() { error("Replace all \"Mat3D\" with \"Mat3\"") } };

struct Mat3DComplex
{ Mat3DComplex() { error("Replace all \"Mat3D\" with \"Mat3\"") } };
