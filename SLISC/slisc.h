// type and class definitions and dependencies
// this header file can be used alone

#pragma once

#ifdef _DEBUG
// this will not check the last index
#define _CHECKBOUNDS_ 1
#endif

// all the system #include's we'll ever need
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include <string.h>

namespace slisc
{

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
typedef const long long int Llong_I; // 64 bit integer
typedef long long int Llong, Llong_O, Llong_IO;
typedef const unsigned long long int Ullong_I;
typedef unsigned long long int Ullong, Ullong_O, Ullong_IO;
#endif

#ifndef _USE_Int_AS_LONG
typedef Llong Long;
#else
typedef Int Long;
#endif
typedef const Long Long_I;
typedef Long Long_O, Long_IO;

typedef const char Char_I; // 8 bit integer
typedef char Char, Char_O, Char_IO;
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

// === constants ===

const Doub PI = 3.14159265358979323;
const Doub E = 2.71828182845904524;
const Comp I(0., 1.);

// report error and pause execution
#define error(str) {std::cout << "error: " << __FILE__ << ": line " << __LINE__ << ": " << str << std::endl; getchar();}

template<class T>
inline void memset(T *dest, const T val, Long_I n)
{
	T *end = dest + n;
	for (; dest < end; ++dest)
		*dest = val;
}

// For cuSLISC project
#ifdef _CUSLISC_
template <class T> class Gvector;
template <class T> class Gmatrix;
template <class T> class Gmat3d;
#endif

// Base Class for vector/matrix
template <class T>
class Vbase
{
protected:
	Long N; // number of elements
	T *p; // pointer to the first element
	inline void move(Vbase &rhs);
public:
	Vbase() : N(0), p(nullptr) {}
	explicit Vbase(Long_I n) : N(n), p(new T[n]) {}
	T* ptr() { return p; } // get pointer
	inline const T* ptr() const { return p; }
	inline Long_I size() const { return N; }
	inline void resize(Long_I n);
	inline T & operator()(Long_I i);
	inline const T & operator()(Long_I i) const;
	inline T& end(); // last element
	inline const T& end() const;
	inline T& end(Long_I i);
	inline const T& end(Long_I i) const;
	~Vbase() { if (p) delete p; }
};

template <class T>
inline void Vbase<T>::resize(Long_I n)
{
	if (n != N) {
		if (p != nullptr) delete[] p;
		N = n;
		p = n > 0 ? new T[n] : nullptr;
	}
}

template <class T>
inline void Vbase<T>::move(Vbase &rhs)
{
	if (p != nullptr) delete[] p;
	N = rhs.N; rhs.N = 0;
	p = rhs.p; rhs.p = nullptr;
}

template <class T>
inline T & Vbase<T>::operator()(Long_I i)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=N)
	error("Vbase subscript out of bounds");
#endif
	return p[i];
}

template <class T>
inline const T & Vbase<T>::operator()(Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=N)
		error("Vbase subscript out of bounds");
#endif
	return p[i];
}

template <class T>
inline T & Vbase<T>::end()
{
#ifdef _CHECKBOUNDS_
	if (N < 1)
		error("Using end() for empty object");
#endif
	return p[N-1];
}

template <class T>
inline const T & Vbase<T>::end() const
{
#ifdef _CHECKBOUNDS_
	if (N < 1)
		error("Using end() for empty object");
#endif
	return p[N-1];
}

template <class T>
inline T& Vbase<T>::end(Long_I i)
{
#ifdef _CHECKBOUNDS_
	if (i <= 0 || i > N)
		error("index out of bound");
#endif
	return p[N-i];
}

template <class T>
inline const T& Vbase<T>::end(Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i <= 0 || i > N)
		error("index out of bound");
#endif
	return p[N-i];
}

// Vector Class

template <class T>
class Vector : public Vbase<T>
{
public:
	typedef Vbase<T> Base;
	using Base::p;
	using Base::N;
	Vector() {}
	explicit Vector(Long_I n): Base(n) {}
	Vector(Long_I n, const T &a) //initialize to constant value
	: Vector(n) { memset(p,a,n); }
	Vector(Long_I n, const T *a) // Initialize to array
	: Vector(n) { memcpy(p, a, n*sizeof(T)); }
	Vector(const Vector &rhs);	// Copy constructor forbidden
	inline Vector & operator=(const Vector &rhs);	// copy assignment
	inline Vector & operator=(const T &rhs);  // assign to constant value
#ifdef _CUSLISC_
	Vector & operator=(const Gvector<T> &rhs) // copy from GPU vector
	{ rhs.get(*this); return *this; }
#endif
	inline void operator<<(Vector &rhs); // move data and rhs.resize(0)
	inline T & operator[](Long_I i);	//i'th element
	inline const T & operator[](Long_I i) const;
	inline void resize(Long_I n) {Base::resize(n);} // resize (contents not preserved)
	template <class T1>
	inline void resize(const Vector<T1> &v) {resize(v.size());}
};

template <class T>
Vector<T>::Vector(const Vector<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference "
		 "argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Vector<T> & Vector<T>::operator=(const Vector<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!");
	resize(rhs);
	memcpy(p, rhs.p, N*sizeof(T));
	return *this;
}

template <class T>
inline Vector<T>& Vector<T>::operator=(const T &rhs)
{
	if (N) memset(p, rhs, N);
	return *this;
}

template <class T>
inline void Vector<T>::operator<<(Vector<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!");
	Base::move(rhs);
}

template <class T>
inline T & Vector<T>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=N)
	error("Vector subscript out of bounds");
#endif
	return p[i];
}

template <class T>
inline const T & Vector<T>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=N)
	error("Vector subscript out of bounds");
#endif
	return p[i];
}

// Matrix Class

template <class T>
class Matrix : public Vbase<T>
{
	typedef Vbase<T> Base;
	using Base::p;
	using Base::N;
private:
	Long nn, mm;
	T **v;
	inline T ** v_alloc();
public:
	Matrix();
	Matrix(Long_I n, Long_I m);
	Matrix(Long_I n, Long_I m, const T &a);	//Initialize to constant
	Matrix(Long_I n, Long_I m, const T *a);	// Initialize to array
	Matrix(const Matrix &rhs);		// Copy constructor
	inline Matrix & operator=(const Matrix &rhs);	// copy assignment
	inline Matrix & operator=(const T &rhs);
#ifdef _CUSLISC_
	Matrix & operator=(const Gmatrix<T> &rhs) // copy from GPU vector
	{ rhs.get(*this); return *this; }
#endif
	inline void operator<<(Matrix &rhs); // move data and rhs.resize(0, 0)
	inline T* operator[](Long_I i);	//subscripting: pointer to row i
	inline const T* operator[](Long_I i) const;
	inline Long nrows() const;
	inline Long ncols() const;
	inline void resize(Long_I n, Long_I m); // resize (contents not preserved)
	template <class T1>
	inline void resize(const Matrix<T1> &a);
	~Matrix();
};

template <class T>
inline T** Matrix<T>::v_alloc()
{
	if (N == 0) return nullptr;
	T **v = new T*[nn];
	v[0] = p;
	for (Long i = 1; i<nn; i++)
		v[i] = v[i-1] + mm;
	return v;
}

template <class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(nullptr) {}

template <class T>
Matrix<T>::Matrix(Long_I n, Long_I m) : Base(n*m), nn(n), mm(m), v(v_alloc()) {}

template <class T>
Matrix<T>::Matrix(Long_I n, Long_I m, const T &s) : Matrix(n, m)
{ memset(p, s, N); }

template <class T>
Matrix<T>::Matrix(Long_I n, Long_I m, const T *ptr) : Matrix(n, m)
{ memcpy(p, ptr, N*sizeof(T)); }

template <class T>
Matrix<T>::Matrix(const Matrix<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!");
	resize(rhs.nn, rhs.mm);
	memcpy(p, rhs.p, N*sizeof(T));
	return *this;
}

template <class T>
inline Matrix<T> & Matrix<T>::operator=(const T &rhs)
{
	if (N) memset(p, rhs, N);
	return *this;
}

template <class T>
inline void Matrix<T>::operator<<(Matrix<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!");
	Base::move(rhs);
	if (v) delete v;
	nn = rhs.nn; mm = rhs.mm; v = rhs.v;
	rhs.nn = rhs.mm = 0; rhs.v = nullptr;;
}

template <class T>
inline T* Matrix<T>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=nn)
		error("Matrix subscript out of bounds");
#endif
	return v[i];
}

template <class T>
inline const T* Matrix<T>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=nn)
		error("Matrix subscript out of bounds");
#endif
	return v[i];
}

template <class T>
inline Long Matrix<T>::nrows() const
{ return nn; }

template <class T>
inline Long Matrix<T>::ncols() const
{ return mm; }

template <class T>
inline void Matrix<T>::resize(Long_I n, Long_I m)
{
	if (n != nn || m != mm) {
		Base::resize(n*m);
		nn = n; mm = m;
		if (v) delete v;
		v = v_alloc();
	}
}

template <class T>
template <class T1>
inline void Matrix<T>::resize(const Matrix<T1> &a)
{ resize(a.nrows(), a.ncols()); }

template <class T>
Matrix<T>::~Matrix()
{ if(v) delete v; }

// 3D Matrix Class

template <class T>
class Mat3d : public Vbase<T>
{
	typedef Vbase<T> Base;
	using Base::p;
	using Base::N;
private:
	Long nn;
	Long mm;
	Long kk;
	T ***v;
	inline T *** v_alloc();
	inline void v_free();
public:
	Mat3d();
	Mat3d(Long_I n, Long_I m, Long_I k);
	Mat3d(Long_I n, Long_I m, Long_I k, const T &a);
	Mat3d(const Mat3d &rhs);   // Copy constructor
	inline Mat3d & operator=(const Mat3d &rhs);	// copy assignment
	inline Mat3d & operator=(const T &rhs);
#ifdef _CUSLISC_
	Mat3d & operator=(const Gmat3d<T> &rhs) // copy from GPU vector
	{ rhs.get(*this); return *this; }
#endif
	inline void operator<<(Mat3d &rhs); // move data and rhs.resize(0, 0, 0)
	inline void resize(Long_I n, Long_I m, Long_I k);
	template <class T1>
	inline void resize(const Mat3d<T1> &a);
	inline T*const * operator[](Long_I i);	//subscripting: pointer to row i
	inline const T*const * operator[](Long_I i) const;
	inline Long dim1() const;
	inline Long dim2() const;
	inline Long dim3() const;
	~Mat3d();
};

template <class T>
inline T *** Mat3d<T>::v_alloc()
{
	if (N == 0) return nullptr;
	Long i;
	Long nnmm = nn*mm;
	T **v0 = new T*[nnmm]; v0[0] = p;
	for (i = 1; i < nnmm; ++i)
		v0[i] = v0[i - 1] + kk;
	T ***v = new T**[nn]; v[0] = v0;
	for(i = 1; i < nn; ++i)
		v[i] = v[i-1] + mm;
	return v;
}

template <class T>
inline void Mat3d<T>::v_free()
{
	if (v != nullptr) {
		delete v[0]; delete v;
	}
}

template <class T>
Mat3d<T>::Mat3d(): nn(0), mm(0), kk(0), v(nullptr) {}

template <class T>
Mat3d<T>::Mat3d(Long_I n, Long_I m, Long_I k) : Base(n*m*k), nn(n), mm(m), kk(k),
	v(v_alloc()) {}

template <class T>
Mat3d<T>::Mat3d(Long_I n, Long_I m, Long_I k, const T &s) : Mat3d(n, m, k)
{ memset(p, s, n*m*k); }

template <class T>
Mat3d<T>::Mat3d(const Mat3d<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Mat3d<T> &Mat3d<T>::operator=(const Mat3d<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!");
	resize(rhs.nn, rhs.mm, rhs.kk);
	memcpy(p, rhs.p, N*sizeof(T));
	return *this;
}

template <class T>
inline Mat3d<T> & Mat3d<T>::operator=(const T &rhs)
{
	if (N) memset(p, rhs, N);
	return *this;
}

template <class T>
inline void Mat3d<T>::operator<<(Mat3d<T> &rhs)
{
	if (this == &rhs) error("self move is forbidden!");
	Base::move(rhs);
	nn = rhs.nn; mm = rhs.mm; kk = rhs.kk;
	v_free(); v = rhs.v;
	rhs.nn = rhs.mm = rhs.kk = 0;
	rhs.v = nullptr;
}

template <class T>
inline void Mat3d<T>::resize(Long_I n, Long_I m, Long_I k)
{
	if (n != nn || m != mm || k != kk) {
		Base::resize(n*m*k);
		nn = n; mm = m; kk = k;
		v_free(); v = v_alloc();
	}
}

template <class T>
template <class T1>
inline void Mat3d<T>::resize(const Mat3d<T1> &a) { resize(a.dim1(), a.dim2(), a.dim3()); }

template <class T>
inline T*const * Mat3d<T>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i >= nn)
		error("Matrix subscript out of bounds");
#endif
	return v[i];
}

template <class T>
inline const T*const * Mat3d<T>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i >= nn)
		error("Matrix subscript out of bounds");
#endif
	return v[i];
}

template <class T>
inline Long Mat3d<T>::dim1() const { return nn; }

template <class T>
inline Long Mat3d<T>::dim2() const { return mm; }

template <class T>
inline Long Mat3d<T>::dim3() const { return kk; }

template <class T>
Mat3d<T>::~Mat3d() { v_free(); }

// Matric and vector types

typedef const Vector<Int> VecInt_I;
typedef Vector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const Vector<Uint> VecUint_I;
typedef Vector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const Vector<Long> VecLong_I;
typedef Vector<Long> VecLong, VecLong_O, VecLong_IO;

typedef const Vector<Llong> VecLlong_I;
typedef Vector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const Vector<Ullong> VecUllong_I;
typedef Vector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const Vector<Char> VecChar_I;
typedef Vector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const Vector<Char*> VecCharp_I;
typedef Vector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const Vector<Uchar> VecUchar_I;
typedef Vector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const Vector<Doub> VecDoub_I;
typedef Vector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const Vector<Doub*> VecDoubp_I;
typedef Vector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const Vector<Comp> VecComp_I;
typedef Vector<Comp> VecComp, VecComp_O, VecComp_IO;

typedef const Vector<Bool> VecBool_I;
typedef Vector<Bool> VecBool, VecBool_O, VecBool_IO;

typedef const Matrix<Int> MatInt_I;
typedef Matrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const Matrix<Uint> MatUint_I;
typedef Matrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const Matrix<Llong> MatLlong_I;
typedef Matrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const Matrix<Ullong> MatUllong_I;
typedef Matrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const Matrix<Char> MatChar_I;
typedef Matrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const Matrix<Uchar> MatUchar_I;
typedef Matrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const Matrix<Doub> MatDoub_I;
typedef Matrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const Matrix<Comp> MatComp_I;
typedef Matrix<Comp> MatComp, MatComp_O, MatComp_IO;

typedef const Matrix<Bool> MatBool_I;
typedef Matrix<Bool> MatBool, MatBool_O, MatBool_IO;

typedef const Mat3d<Doub> Mat3Doub_I;
typedef Mat3d<Doub> Mat3Doub, Mat3Doub_O, Mat3Doub_IO;

typedef const Mat3d<Comp> Mat3Comp_I;
typedef Mat3d<Comp> Mat3Comp, Mat3Comp_O, Mat3Comp_IO;

// macro-like functions (don't use them in your code ever, write similar utilities in "algorithm.h")

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

} // namespace slisc
