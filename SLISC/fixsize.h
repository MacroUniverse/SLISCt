// fixed-size vector/matrix classes

#include "arithmetic.h"

namespace slisc {

template <class T, Long N> class FixVec;
template <class T, Long Nr, Long Nc> class FixCmat;

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

// Base Class for vector/matrix
template <class T, Long N>
class FixVec
{
protected:
	T m_data[N]; // pointer to the first element
public:
	typedef T value_type;
	FixVec() {}
	FixVec(const T &s) { *this = s; }
	T* ptr() { return m_data; } // get pointer
	const T* ptr() const { return m_data; }
	Long size() const { return N; }
	T & operator[](Long_I i);
	const T & operator[](Long_I i) const;
	T & operator()(Long_I i);
	const T & operator()(Long_I i) const;
	T& end(Long_I i = 1);
	const T& end(Long_I i = 1) const;
	FixVec & operator=(const FixVec &rhs);
	template <class T1>
	FixVec & operator=(const FixVec<T1, N> &rhs);
	FixVec & operator=(const T &rhs); // for scalar
};

template <class T, Long N>
inline T & FixVec<T, N>::operator[](Long_I i)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=N)
		error("FixVec subscript out of bounds");
#endif
	return m_data[i];
}

template <class T, Long N>
inline const T & FixVec<T, N>::operator[](Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=N)
		error("FixVec subscript out of bounds");
#endif
	return m_data[i];
}

template <class T, Long N>
inline T & FixVec<T, N>::operator()(Long_I i)
{ return (*this)[i]; }

template <class T, Long N>
inline const T & FixVec<T, N>::operator()(Long_I i) const
{ return (*this)[i]; }

template <class T, Long N>
inline FixVec<T, N> & FixVec<T, N>::operator=(const FixVec<T, N> &rhs)
{ return operator=<T>(rhs); }

template <class T, Long N>
inline FixVec<T, N> & FixVec<T, N>::operator=(const T &rhs)
{
	vecset(m_data, rhs, N);
	return *this;
}

template <class T, Long N> template <class T1>
inline FixVec<T, N> & FixVec<T, N>::operator=(const FixVec<T1, N> &rhs)
{
	veccpy(m_data, rhs.ptr(), N);
	return *this;
}

template <class T, Long N>
inline T & FixVec<T, N>::end(Long_I i)
{
#ifdef _CHECKBOUNDS_
	if (i <= 0 || i > N)
		error("index out of bound");
#endif
	return m_data[N-i];
}

template <class T, Long N>
inline const T & FixVec<T, N>::end(Long_I i) const
{
#ifdef _CHECKBOUNDS_
	if (i <= 0 || i > N)
		error("index out of bound");
#endif
	return m_data[N-i];
}

// Column major Matrix Class

template <class T, Long Nr, Long Nc>
class FixCmat : public FixVec<T, Nr*Nc>
{
private:
	typedef FixVec<T, Nr*Nc> Base;
	using Base::m_data;
public:
	using Base::operator();
	using Base::operator=;
	FixCmat() {}
	FixCmat(const T &s); //Initialize to constant
	FixCmat(const FixCmat &rhs); // Copy constructor
	FixCmat & operator=(const T &rhs);
	FixCmat & operator=(const FixCmat &rhs); // copy assignment
	template <class T1, Long Nr1, Long Nc1>
	FixCmat & operator=(const FixCmat<T1, Nr1, Nc1> &rhs);
	T& operator()(Long_I i, Long_I j);	// double indexing
	const T& operator()(Long_I i, Long_I j) const;
	Long nrows() const;
	Long ncols() const;
};

template <class T, Long Nr, Long Nc>
FixCmat<T, Nr, Nc>::FixCmat(const T &s)
{ *this = s; }

template <class T, Long Nr, Long Nc>
FixCmat<T, Nr, Nc>::FixCmat(const FixCmat<T, Nr, Nc> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T, Long Nr, Long Nc>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const FixCmat<T, Nr, Nc> &rhs)
{
	return operator=<T, Nr, Nc>(rhs);
}

template <class T, Long Nr, Long Nc>
template <class T1, Long Nr1, Long Nc1>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const FixCmat<T1, Nr1, Nc1> &rhs)
{
	Nr = rhs.nrows();
	Nc = rhs.ncols();
	Base::operator=(rhs);
	return *this;
}

template <class T, Long Nr, Long Nc>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const T &rhs)
{
	Base::operator=(rhs);
	return *this;
}

template <class T, Long Nr, Long Nc>
inline T & FixCmat<T, Nr, Nc>::operator()(Long_I i, Long_I j)
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= Nr || j < 0 || j >= Nc)
		error("Matrix subscript out of bounds");
#endif
	return m_data[i+Nr*j];
}

template <class T, Long Nr, Long Nc>
inline const T & FixCmat<T, Nr, Nc>::operator()(Long_I i, Long_I j) const
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= Nr || j < 0 || j >= Nc)
		error("Matrix subscript out of bounds");
#endif
	return m_data[i+Nr*j];
}

template <class T, Long Nr, Long Nc>
inline Long FixCmat<T, Nr, Nc>::nrows() const
{ return Nr; }

template <class T, Long Nr, Long Nc>
inline Long FixCmat<T, Nr, Nc>::ncols() const
{ return Nc; }

// arithmetics

template <class T, Long N, class T1, Long N1>
Bool operator==(const FixVec<T, N> &v, const FixVec<T1, N1> &v1)
{ equals_to_vv(v, v1); }

template <class T, Long N, class T1, Long N1>
Bool operator!=(const FixVec<T, N> &v, const FixVec<T1, N1> &v1)
{ return !(v==v1); }

template <class T, Long N, class T1>
Bool operator==(const FixVec<T, N> &v, const T1 &s)
{ return equals_to_vs(v.ptr(), s, v.size()); }

template <class T, Long N, class T1>
Bool operator!=(const FixVec<T, N> &v, const T1 &s)
{ return !(v==s); }

template <class T, Long Nr, Long Nc, class T1, Long Nr1, Long Nc1>
Bool operator==(const FixCmat<T, Nr, Nc> &v, const FixCmat<T1, Nr, Nc> &v1)
{ return equals_to_vv(v, v1); }

template <class T, Long Nr, Long Nc, class T1, Long Nr1, Long Nc1>
Bool operator!=(const FixCmat<T, Nr, Nc> &v, const FixCmat<T1, Nr, Nc> &v1)
{ return !(v==v1); }

template <class T, Long Nr, Long Nc, class T1>
Bool operator==(const FixCmat<T, Nr, Nc> &v, const T1 &s)
{ return equals_to_vs(v.ptr(), s, v.size()); }

template <class T, Long Nr, Long Nc, class T1>
Bool operator!=(const FixCmat<T, Nr, Nc> &v, const T1 &s)
{ return !(v==s); }

}
