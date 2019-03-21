#pragma once
#include "vector.h"

namespace slisc {

// contiguous slice vector class
template <class T>
class Svector : public Vector<T>
{
public:
	typedef Vector<T> Base;
	using Base::m_p;
	using Base::m_N;
	Svector();
	Svector(Long_I N);
	Svector(const T *ptr, Long_I N);
	// There is no upper bound checking of N, use with care
	void set_N(Long_I N);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I N);
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
inline void Svector<T>::set_N(Long_I N)
{
#ifdef SLS_CHECK_SHAPE
	if (N <= 0) error("illegal N!");
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
inline Svector<T>::~Svector()
{
	m_p = nullptr; m_N = 0;
}

}
