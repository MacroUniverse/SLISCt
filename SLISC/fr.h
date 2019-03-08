// 3d function in spherical harmonic expansion
#pragma once
#include "scalar_arith.h"
#include "cmat.h"

namespace slisc {

template <class T>
class FR : public Cmat<Vector<T>>
{
private:
	Long m_Mmin, m_Mmax;
public:
	typedef Cmat<Vector<T>> Base;
	using Base::nrows;
	using Base::ncols;
	using Base::operator();
	FR(Long_I Lmax, Long_I Nr); // L = 0...NL-1, M = 0
	FR(Long_I Lmax, Long_I Mmin, Long_I Mmax, Long_I Nr); // L = 0...NL-1, M = Mmin...Mmax
	Long Lmax() const;
	Long Mmin() const;
	Long Mmax() const;
	// check if fr{L,M} exists, for any L, M
	Bool exist(Long_I L, Long_I M) const;
	const Vector<T> &get(Long_I L, Long_I M) const;
	Vector<T> &get(Long_I L, Long_I M);
	//const T get(Long_I L, Long_I M, Long_I i);
	~FR();
};

template <class T>
FR<T>::FR(Long_I Lmax, Long_I Nr) :
	Base(Lmax+1,1), m_Mmin(0), m_Mmax(0)
{
	for (Long i = 0; i <= Lmax; ++i)
		(*this)[i].resize(Nr);
}

template <class T>
FR<T>::FR(Long_I Lmax, Long_I Mmin, Long_I Mmax, Long_I Nr) :
	Base(Lmax + 1, Mmax - Mmin + 1), m_Mmin(Mmin), m_Mmax(Mmax)
{
	for (Long L = 0; L <= Lmax; ++L)
		for (Long M = max(-L, Mmin); M <= min(L,Mmax); ++M)
			(*this)(L, M-Mmin).resize(Nr);
}

template <class T>
Long FR<T>::Lmax() const
{
	return nrows() - 1;
}

template <class T>
Long FR<T>::Mmin() const
{
	return m_Mmin;
}

template <class T>
Long FR<T>::Mmax() const
{
	return m_Mmax;
}

template <class T>
Bool FR<T>::exist(Long_I L, Long_I M) const
{
#ifdef SLS_CHECK_BOUNDS
	if (L < 0) error("L < 0 is illegal!");
#endif
	if (L <= Lmax() && M >= m_Mmin && M <= m_Mmax && get(L, M).size() > 0) return true;
	return false;
}

template <class T>
const Vector<T> &FR<T>::get(Long_I L, Long_I M) const
{
#ifdef SLS_CHECK_BOUNDS
	if (L < 0 || L > Lmax() || M < m_Mmin || M > m_Mmax)
		error("{L,M} out of bound!");
#endif
	return (*this)(L, M - m_Mmin);
}

template <class T>
Vector<T> &FR<T>::get(Long_I L, Long_I M)
{
#ifdef SLS_CHECK_BOUNDS
	if (!exist(L, M)) error("{L,M} out of bound!");
#endif
	return (*this)(L, M - m_Mmin);
}

template <class T>
FR<T>::~FR()
{
	for (Long i = 0; i < nrows(); ++i)
		for (Long j = 0; j < ncols(); ++j)
			(*this)(i, j).resize(0);
}
}
