// vector container
#pragma once
#include "meta.h"

#ifndef NDEBUG
#define SLS_CHECK_BOUNDS
#define SLS_CHECK_SHAPE
#endif

namespace slisc {

// array copying
template<class T>
inline void vecset(T *dest, const T &val, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = val;
}

template<class T>
inline void veccpy(T *dest, const T *src, Long_I n)
{
	memcpy(dest, src, n * sizeof(T));
}

template<class T, class T1, SLS_IF(is_promo<T,T1>())>
inline void veccpy(T *dest, const T1 *src, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = src[i];
}

// Base Class for vector/matrix
template <class T>
class Vbase
{
protected:
	Long m_N; // number of elements
	T *m_p; // pointer to the first element
public:
	typedef T value_type;
	Vbase() : m_N(0), m_p(nullptr) {}
	explicit Vbase(Long_I N) : m_N(N), m_p(new T[N]) {}
	T* ptr() { return m_p; } // get pointer
	const T* ptr() const { return m_p; }
	Long size() const { return m_N; }
	void resize(Long_I N);
	T & operator[](Long_I i);
	const T & operator[](Long_I i) const;
	T & operator()(Long_I i);
	const T & operator()(Long_I i) const;
	T& end(Long_I i = 1);
	const T& end(Long_I i = 1) const;
	Vbase & operator=(const Vbase &rhs);
	template <class T1>
	Vbase & operator=(const Vbase<T1> &rhs);
	Vbase & operator=(const T &rhs); // for scalar
	void operator<<(Vbase &rhs);
	~Vbase() {
		if (m_p)
			delete m_p;
	}
};

template <class T>
inline void Vbase<T>::resize(Long_I N)
{
	if (N != m_N) {
		if (m_p != nullptr)
			delete[] m_p;
		m_N = N;
		m_p = N > 0 ? new T[N] : nullptr;
	}
}

template <class T>
inline void Vbase<T>::operator<<(Vbase &rhs)
{
	if (this == &rhs)
		error("self move is forbidden!");
	if (m_p != nullptr) delete[] m_p;
	m_N = rhs.m_N; rhs.m_N = 0;
	m_p = rhs.m_p; rhs.m_p = nullptr;
}

template <class T>
inline T & Vbase<T>::operator[](Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
if (i<0 || i>=m_N)
	error("Vbase subscript out of bounds");
#endif
	return m_p[i];
}

template <class T>
inline const T & Vbase<T>::operator[](Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i<0 || i>=m_N)
		error("Vbase subscript out of bounds");
#endif
	return m_p[i];
}

template <class T>
inline T & Vbase<T>::operator()(Long_I i)
{ return (*this)[i]; }

template <class T>
inline const T & Vbase<T>::operator()(Long_I i) const
{ return (*this)[i]; }

template <class T>
inline Vbase<T> & Vbase<T>::operator=(const Vbase<T> &rhs)
{ return operator=<T>(rhs); }

template <class T>
inline Vbase<T> & Vbase<T>::operator=(const T &rhs)
{
	vecset(m_p, rhs, m_N);
	return *this;
}

template <class T> template <class T1>
inline Vbase<T> & Vbase<T>::operator=(const Vbase<T1> &rhs)
{
	resize(rhs.size());
	veccpy(m_p, rhs.ptr(), m_N);
	return *this;
}

template <class T>
inline T & Vbase<T>::end(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
	if (i <= 0 || i > m_N)
		error("index out of bound");
#endif
	return m_p[m_N-i];
}

template <class T>
inline const T & Vbase<T>::end(Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i <= 0 || i > m_N)
		error("index out of bound");
#endif
	return m_p[m_N-i];
}

// Vector Class

template <class T>
class Vector : public Vbase<T>
{
protected:
	typedef Vbase<T> Base;
	using Base::m_p;
	using Base::m_N;
public:
	using Base::resize;
	using Base::operator=;
	Vector() {}
	explicit Vector(Long_I N): Base(N) {}
	Vector(Long_I N, const T &s) //initialize to constant value
	: Vector(N) { *this = s; }
	Vector(Long_I N, const T *a) // Initialize to array
	: Vector(N) { veccpy(m_p, a, N); }
	Vector(const Vector &rhs);	// Copy constructor forbidden
	static constexpr Int ndims() { return 1; }
	Vector &operator=(const Vector &rhs);
	template <class T1>
	Vector &operator=(const Vector<T1> &rhs);
#ifdef _CUSLISC_
	Vector & operator=(const Gvector<T> &rhs) // copy from GPU vector
	{ rhs.get(*this); return *this; }
#endif
	void operator<<(Vector &rhs); // move data and rhs.resize(0)
	template <class T1>
	void resize(const Vector<T1> &v) {resize(v.size());}
};

template <class T>
Vector<T>::Vector(const Vector<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference "
		 "argument for function input or output, and use \"=\" to copy!");
}

template <class T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs)
{
	return operator=<T>(rhs);
}

template <class T> template <class T1>
Vector<T> &Vector<T>::operator=(const Vector<T1> &rhs)
{
	Base::operator=(rhs);
	return *this;
}

template <class T>
inline void Vector<T>::operator<<(Vector<T> &rhs)
{
	Base::operator<<(rhs);
}

} // namespace slisc
