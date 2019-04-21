#pragma once
#include "vector.h"

namespace slisc {

// contiguous slice vector class
template <class T>
class Svector : public Vector<T>
{
public:
	typedef Vector<T> Base;
	using Base::value_type;
	using Base::m_p;
	using Base::m_N;
	Svector();
	Svector(Long_I N);
	Svector(const T *ptr, Long_I N);
	// There is no upper bound checking of N, use with care
	void set_size(Long_I N);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I N);
	void resize(Long_I N);
	~Svector();
};

template <class T>
inline Svector<T>::Svector()
{}

template <class T>
inline Svector<T>::Svector(Long_I N)
{
	m_N = N;
}

template <class T>
inline Svector<T>::Svector(const T *ptr, Long_I N)
{
	// TODO: might be inefficient since m_p, m_N are initialized to 0 by Vector class first
	// a constructor of Vbase/Vector that leaves things uninitialized might be added.
	// same problem with other constructors and destructor
	m_p = (T *)ptr; m_N = N;
}

template<class T>
inline void Svector<T>::set_size(Long_I N)
{
#ifdef SLS_CHECK_SHAPE
	if (N <= 0) SLS_ERR("illegal N!");
#endif
	m_N = N;
}

template<class T>
inline void Svector<T>::set_ptr(const T * ptr)
{
	m_p = (T *)ptr;
}

template<class T>
inline void Svector<T>::set(const T * ptr, Long_I N)
{
	m_p = (T *)ptr; m_N = N;
}

template<class T>
inline void Svector<T>::resize(Long_I N)
{
	SLS_ERR("Svector cannot be resized, consider using set_size()!");
}

template<class T>
inline Svector<T>::~Svector()
{
	m_p = nullptr; m_N = 0;
}

// contiguous slice matrix class (column major)
template <class T>
class Scmat : public Cmat<T>
{
public:
	typedef Cmat<T> Base;
	using Base::value_type;
	using Base::m_p;
	using Base::m_N;
	using Base::m_Nr;
	using Base::m_Nc;
	Scmat();
	Scmat(Long_I Nr, Long_I Nc);
	Scmat(const T *ptr, Long_I Nr, Long_I Nc);
	// There is no upper bound checking of N, use with care
	void set_size(Long_I Nr, Long_I Nc);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I Nr, Long_I Nc);
	void resize(Long_I Nr, Long_I Nc); // forbidden
	~Scmat();
};

template <class T>
inline Scmat<T>::Scmat()
{}

template <class T>
inline Scmat<T>::Scmat(Long_I Nr, Long_I Nc)
{
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline Scmat<T>::Scmat(const T *ptr, Long_I Nr, Long_I Nc)
{
	// TODO: might be inefficient since m_p, m_N are initialized to 0 by Vector class first
	// a constructor of Vbase/Vector that leaves things uninitialized might be added.
	// same problem with other constructors and destructor
	m_p = (T *)ptr;
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template<class T>
inline void Scmat<T>::set_size(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr <= 0 || Nc <= 0) SLS_ERR("illegal Nr or Nc!");
#endif
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template<class T>
inline void Scmat<T>::set_ptr(const T * ptr)
{
	m_p = (T *)ptr;
}

template<class T>
inline void Scmat<T>::set(const T * ptr, Long_I Nr, Long_I Nc)
{
	m_p = (T *)ptr;
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template<class T>
inline void Scmat<T>::resize(Long_I Nr, Long_I Nc)
{
	SLS_ERR("Scmat<> cannot be resized, consider using set_size()!");
}

template<class T>
inline Scmat<T>::~Scmat()
{
	m_p = nullptr;
}

}
