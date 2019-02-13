// column-major matrix container
#pragma once
#include "matrix.h"

namespace slisc {

// Column major Matrix Class
// don't derive from Matrix<T>

template <class T>
class Cmat : public Vbase<T>
{
private:
	typedef Vbase<T> Base;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc;
public:
	using Base::operator();
	using Base::operator=;
	Cmat();
	Cmat(Long_I Nr, Long_I Nc);
	Cmat(Long_I Nr, Long_I Nc, const T &s);	//Initialize to constant
	Cmat(Long_I Nr, Long_I Nc, const T *ptr);	// Initialize to array
	Cmat(const Cmat &rhs);		// Copy constructor
	static constexpr Int ndims() { return 2; } // matrix is 2 dimensional
	static constexpr Char major() { return 'c'; } // row major memory
	Cmat & operator=(const Cmat &rhs);	// copy assignment
	template <class T1>
	Cmat & operator=(const Cmat<T1> &rhs);
	Cmat & operator=(const T &rhs);
	template <class T1>
	Cmat & operator=(const MatCoo<T1> &rhs);
	template <class T1>
	Cmat & operator=(const MatCooH<T1> &rhs);
	void operator<<(Cmat &rhs); // move data and rhs.resize(0, 0)
	T& operator()(Long_I i, Long_I j);	// double indexing
	const T& operator()(Long_I i, Long_I j) const;
	Long nrows() const;
	Long ncols() const;
	void resize(Long_I Nr, Long_I Nc); // resize (contents not preserved)
	template <class T1>
	void resize(const Cmat<T1> &a);
};

template <class T>
Cmat<T>::Cmat() : m_Nr(0), m_Nc(0) {}

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc) : Base(Nr*Nc), m_Nr(Nr), m_Nc(Nc) {}

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc, const T &s) : Cmat(Nr, Nc)
{ *this = s; }

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc, const T *ptr) : Cmat(Nr, Nc)
{ memcpy(m_p, ptr, m_N*sizeof(T)); }

template <class T>
Cmat<T>::Cmat(const Cmat<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Cmat<T> & Cmat<T>::operator=(const Cmat<T> &rhs)
{
	return operator=<T>(rhs);
}

template <class T> template <class T1>
inline Cmat<T> & Cmat<T>::operator=(const Cmat<T1> &rhs)
{
	m_Nr = rhs.nrows();
	m_Nc = rhs.ncols();
	Base::operator=(rhs);
	return *this;
}

template <class T>
inline Cmat<T> & Cmat<T>::operator=(const T &rhs)
{
	Base::operator=(rhs);
	return *this;
}

template <class T> template <class T1>
inline Cmat<T> & Cmat<T>::operator=(const MatCoo<T1> &rhs)
{
	return coo2dense(*this, rhs);
}

template <class T> template <class T1>
inline Cmat<T> & Cmat<T>::operator=(const MatCooH<T1> &rhs)
{
	return cooh2dense(*this, rhs);
}

template <class T>
inline void Cmat<T>::operator<<(Cmat<T> &rhs)
{
	m_Nr = rhs.m_Nr; m_Nc = rhs.m_Nc;
	rhs.m_Nr = rhs.m_Nc = 0;
	Base::operator<<(rhs);
}

template <class T>
inline T & Cmat<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		error("Matrix subscript out of bounds");
#endif
	return m_p[i+m_Nr*j];
}

template <class T>
inline const T & Cmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		error("Matrix subscript out of bounds");
#endif
	return m_p[i+m_Nr*j];
}

template <class T>
inline Long Cmat<T>::nrows() const
{ return m_Nr; }

template <class T>
inline Long Cmat<T>::ncols() const
{ return m_Nc; }

template <class T>
inline void Cmat<T>::resize(Long_I Nr, Long_I Nc)
{
	if (Nr != m_Nr || Nc != m_Nc) {
		Base::resize(Nr*Nc);
		m_Nr = Nr; m_Nc = Nc;
	}
}

template <class T>
template <class T1>
inline void Cmat<T>::resize(const Cmat<T1> &a)
{ resize(a.nrows(), a.ncols()); }

}
