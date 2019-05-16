#pragma once
#include "vector.h"

namespace slisc {

// contiguous slice vector class
template <class T>
class Svector
{
protected:
	Vector<T> &cast();
	const Vector<T> &cast() const;
public:
	typedef T value_type;
	T *m_p;
	Long m_N;
	Svector();
	Svector(Long_I N);
	Svector(const T *ptr, Long_I N);

	// === Vbase<T> functions ===
	T* ptr(); // get pointer
	const T* ptr() const;
	Long size() const;
	// resize() is a bad idea, don't try to create it!
	T & operator[](Long_I i);
	const T & operator[](Long_I i) const;
	T & operator()(Long_I i);
	const T & operator()(Long_I i) const;
	T& end();
	const T& end() const;
	T& end(Long_I i);
	const T& end(Long_I i) const;
	Svector & operator=(const Svector &rhs);
	template <class T1>
	Svector & operator=(const Svector<T1> &rhs);
	Svector & operator=(const T &rhs); // for scalar

	// === other member functions ===
	// There is no upper bound checking for set_size(), use with care
	void set_size(Long_I N);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I N);
	
	~Svector();
};

template<class T>
inline Vector<T>& Svector<T>::cast()
{
	return *reinterpret_cast<Vector<T>*>(this);
}

template<class T>
inline const Vector<T>& Svector<T>::cast() const
{
	return *reinterpret_cast<const Vector<T>*>(this);
}

template <class T>
inline Svector<T>::Svector() {}

template <class T>
inline Svector<T>::Svector(Long_I N) : m_N(N) {}

template <class T>
inline Svector<T>::Svector(const T *ptr, Long_I N)
	: m_p((T *)ptr), m_N(N) {}

template<class T>
inline T * Svector<T>::ptr()
{
	return cast().ptr();
}

template<class T>
inline const T * Svector<T>::ptr() const
{
	return cast().ptr();
}

template<class T>
inline Long Svector<T>::size() const
{
	return m_N;
}

template<class T>
inline T & Svector<T>::operator[](Long_I i)
{
	return cast()[i];
}

template<class T>
inline const T & Svector<T>::operator[](Long_I i) const
{
	return cast()[i];
}

template<class T>
inline T & Svector<T>::operator()(Long_I i)
{
	return cast()(i);
}

template<class T>
inline const T & Svector<T>::operator()(Long_I i) const
{
	return cast()(i);
}

template<class T>
inline T & Svector<T>::end()
{
	return cast().end();
}

template<class T>
inline const T & Svector<T>::end() const
{
	return cast().end();
}

template<class T>
inline T & Svector<T>::end(Long_I i)
{
	return cast().end(i);
}

template<class T>
inline const T & Svector<T>::end(Long_I i) const
{
	return cast().end(i);
}

template <class T>
inline Svector<T> & Svector<T>::operator=(const Svector<T> &rhs)
{
	return cast() = rhs;
}

template <class T>
inline Svector<T> & Svector<T>::operator=(const T &rhs)
{
	return cast() = rhs;
}

template <class T> template <class T1>
inline Svector<T> & Svector<T>::operator=(const Svector<T1> &rhs)
{
	return cast() = rhs;
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
inline Svector<T>::~Svector()
{}

// contiguous slice matrix class (column major)
template <class T>
class Scmat : public Svector<T>
{
private:
	Cmat<T> &cast();
	const Cmat<T> &cast() const;
public:
	typedef Svector<T> Base;
	using Base::value_type;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc;
	Scmat();
	Scmat(Long_I Nr, Long_I Nc, Long_I N = Nr*Nc);
	Scmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I N = Nr * Nc);

	// === Cmat functions ===
	static constexpr Int ndims() { return 2; } // matrix is 2 dimensional
	static constexpr Char major() { return 'c'; } // row major memory
	Scmat & operator=(const Scmat &rhs);	// copy assignment
	template <class T1>
	Scmat & operator=(const Scmat<T1> &rhs);
	Scmat & operator=(const T &rhs);
	template <class T1>
	Scmat & operator=(const MatCoo<T1> &rhs);
	template <class T1>
	Scmat & operator=(const MatCooH<T1> &rhs);
	T& operator()(Long_I i, Long_I j);	// double indexing
	const T& operator()(Long_I i, Long_I j) const;
	Long nrows() const;
	Long ncols() const;

	// resize() is a bad idea, don't try to create it!

	// There is no upper bound checking of N, use with care
	void set_size(Long_I Nr, Long_I Nc);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I Nr, Long_I Nc);
	~Scmat();
};

template<class T>
inline Cmat<T>& Scmat<T>::cast()
{
	return *reinterpret_cast<Cmat<T>*>(this);
}

template<class T>
inline const Cmat<T>& Scmat<T>::cast() const
{
	return *reinterpret_cast<const Cmat<T>*>(this);
}

template <class T>
inline Scmat<T>::Scmat()
{}

template <class T>
inline Scmat<T>::Scmat(Long_I Nr, Long_I Nc, Long_I N)
	: m_Nr(Nr), m_Nc(Nc), Base(N)
{}

template <class T>
inline Scmat<T>::Scmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I N)
	: Scmat(Nr, Nc, N)
{
	m_p = (T *)ptr;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const Scmat<T> &rhs)
{
	return cast() = rhs;
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const Scmat<T1> &rhs)
{
	return cast() = rhs;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const T &rhs)
{
	return cast() = rhs;
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCoo<T1> &rhs)
{
	return cast() = rhs;
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCooH<T1> &rhs)
{
	return cast() = rhs;
}

template<class T>
inline T & Scmat<T>::operator()(Long_I i, Long_I j)
{
	return cast()(i, j);
}

template<class T>
inline const T & Scmat<T>::operator()(Long_I i, Long_I j) const
{
	return cast()(i, j);
}

template<class T>
inline Long Scmat<T>::nrows() const
{
	return m_Nr;
}

template<class T>
inline Long Scmat<T>::ncols() const
{
	return m_Nc;
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
inline Scmat<T>::~Scmat() {}

} // namespace slisc
