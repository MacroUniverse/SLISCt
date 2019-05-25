#pragma once
#include "vector.h"

namespace slisc {

// contiguous slice vector class
template <class T>
class Svector
{
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
	template <class Tv, SLS_IF(is_dense_vec<Tv>())>
	Svector & operator=(const Tv &rhs);
	Svector & operator=(const T &rhs); // for scalar

	// === other member functions ===
	// There is no bound checking, use with care
	void set_size(Long_I N);
	void set_ptr(const T *ptr);
	void set(const T *ptr, Long_I N);
	void next(); // m_ptr += m_N
	void last(); // m_ptr -= m_N
	void shift(Long_I N); // m_ptr += N;
	
	~Svector();
};

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
#ifdef SLS_CHECK_BOUNDS
	if (m_N == 0)
		SLS_ERR("using ptr() for empty container!");
#endif
	return m_p;
}

template<class T>
inline const T * Svector<T>::ptr() const
{
#ifdef SLS_CHECK_BOUNDS
	if (m_N == 0)
		SLS_ERR("using ptr() for empty container!");
#endif
	return m_p;
}

template<class T>
inline Long Svector<T>::size() const
{
	return m_N;
}

template<class T>
inline T & Svector<T>::operator[](Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_N)
		SLS_ERR("Vbase subscript out of bounds");
#endif
	return m_p[i];
}

template<class T>
inline const T & Svector<T>::operator[](Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_N)
		SLS_ERR("Vbase subscript out of bounds");
#endif
	return m_p[i];
}

template<class T>
inline T & Svector<T>::operator()(Long_I i)
{
	return (*this)[i];
}

template<class T>
inline const T & Svector<T>::operator()(Long_I i) const
{
	return (*this)[i];
}

template<class T>
inline T & Svector<T>::end()
{
#ifdef SLS_CHECK_BOUNDS
	if (m_N == 0)
		SLS_ERR("tring to use end() on empty vector!");
#endif
	return m_p[m_N - 1];
}

template<class T>
inline const T & Svector<T>::end() const
{
#ifdef SLS_CHECK_BOUNDS
	if (m_N == 0)
		SLS_ERR("tring to use end() on empty vector!");
#endif
	return m_p[m_N - 1];
}

template<class T>
inline T & Svector<T>::end(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
	if (i <= 0 || i > m_N)
		SLS_ERR("index out of bound");
#endif
	return m_p[m_N - i];
}

template<class T>
inline const T & Svector<T>::end(Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i <= 0 || i > m_N)
		SLS_ERR("index out of bound");
#endif
	return m_p[m_N - i];
}

template <class T>
inline Svector<T> & Svector<T>::operator=(const Svector<T> &rhs)
{
	veccpy(m_p, rhs.ptr(), m_N);
	return *this;
}

template <class T> template <class Tv, SLS_IF0(is_dense_vec<Tv>())>
inline Svector<T> & Svector<T>::operator=(const Tv &rhs)
{
	veccpy(m_p, rhs.ptr(), m_N);
	return *this;
}

template <class T>
inline Svector<T> & Svector<T>::operator=(const T &rhs)
{
	vecset(m_p, rhs, m_N);
	return *this;
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
inline void Svector<T>::next()
{
	m_p += m_N;
}

template<class T>
inline void Svector<T>::last()
{
	m_p -= m_N;
}

template<class T>
inline void Svector<T>::shift(Long_I N)
{
	m_p += N;
}

template<class T>
inline Svector<T>::~Svector()
{}

// contiguous slice matrix class (column major)
template <class T>
class Scmat : public Svector<T>
{
public:
	typedef Svector<T> Base;
	using Base::value_type;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc;
	Scmat();
	Scmat(Long_I Nr, Long_I Nc);
	Scmat(const T *ptr, Long_I Nr, Long_I Nc);

	// === Cmat functions ===
	Scmat & operator=(const Scmat &rhs);	// copy assignment
	template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
	Scmat & operator=(const Tmat &rhs);
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
	void next(); // m_ptr += m_N
	void last(); // m_ptr -= m_N
	void shift(Long_I N); // m_ptr += N;
	~Scmat();
};

template <class T>
inline Scmat<T>::Scmat() {}

template <class T>
inline Scmat<T>::Scmat(Long_I Nr, Long_I Nc)
	: m_Nr(Nr), m_Nc(Nc), Base(Nr*Nc) {}

template <class T>
inline Scmat<T>::Scmat(const T *ptr, Long_I Nr, Long_I Nc)
	: Scmat(Nr, Nc)
{
	m_p = (T *)ptr;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const Scmat<T> &rhs)
{
	copy(*this, rhs);
	return *this;
}

template <class T> template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
inline Scmat<T> & Scmat<T>::operator=(const Tmat &rhs)
{
	copy(*this, rhs);
	return *this;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const T &rhs)
{
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCoo<T1> &rhs)
{
	return coo2dense(*this, rhs);
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCooH<T1> &rhs)
{
	return cooh2dense(*this, rhs);
}

template <class T>
inline T & Scmat<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("Matrix subscript out of bounds");
#endif
	return m_p[i + m_Nr * j];
}

template <class T>
inline const T & Scmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("Matrix subscript out of bounds");
#endif
	return m_p[i + m_Nr * j];
}

template <class T>
inline Long Scmat<T>::nrows() const
{
	return m_Nr;
}

template <class T>
inline Long Scmat<T>::ncols() const
{
	return m_Nc;
}

template <class T>
inline void Scmat<T>::set_size(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr <= 0 || Nc <= 0) SLS_ERR("illegal Nr or Nc!");
#endif
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat<T>::set_ptr(const T * ptr)
{
	m_p = (T *)ptr;
}

template <class T>
inline void Scmat<T>::set(const T * ptr, Long_I Nr, Long_I Nc)
{
	m_p = (T *)ptr;
	m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat<T>::next()
{
	m_p += m_N;
}

template <class T>
inline void Scmat<T>::last()
{
	m_p -= m_N;
}

template <class T>
inline void Scmat<T>::shift(Long_I N)
{
	m_p += N;
}

template <class T>
inline Scmat<T>::~Scmat() {}

template <class T>
class Dcmat
{
private:
	T *m_p;
	Long m_N;
	Long m_Nr, m_Nc;
	Long m_lda;
public:
	Dcmat();
	Dcmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I lda);
	void set(const T *ptr, Long_I Nr, Long_I Nc, Long_I lda);

	// === Cmat member functions ===
	//Cmat & operator=(const Cmat &rhs);	// copy assignment
	//template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
	//Cmat & operator=(const Tmat &rhs);
	//Cmat & operator=(const T &rhs);
	//template <class T1>
	//Cmat & operator=(const MatCoo<T1> &rhs);
	//template <class T1>
	//Cmat & operator=(const MatCooH<T1> &rhs);
	T& operator()(Long_I i, Long_I j);	// double indexing
	const T& operator()(Long_I i, Long_I j) const;
	Long nrows() const;
	Long ncols() const;
	Long lda() const;
};

template <class T>
Dcmat<T>::Dcmat() {}

template <class T>
Dcmat<T>::Dcmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I lda)
	: m_p(ptr), m_Nr(Nr), m_Nc(Nc), m_lda(lda)
{}

template <class T>
void Dcmat<T>::set(const T *ptr, Long_I Nr, Long_I Nc, Long_I lda)
{
	m_p = ptr; m_Nr = Nr; m_Nc = Nc; m_lda = lda;
}

template <class T>
inline T & Dcmat<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("Matrix subscript out of bounds");
#endif
	return m_p[i + m_lda * j];
}

template <class T>
inline const T & Dcmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("Matrix subscript out of bounds");
#endif
	return m_p[i + m_lda * j];
}

template <class T>
inline Long Dcmat<T>::nrows() const
{
	return m_Nr;
}

template <class T>
inline Long Dcmat<T>::ncols() const
{
	return m_Nc;
}

template <class T>
inline Long Dcmat<T>::lda() const
{
	return m_lda;
}

// === arithmetics ===

// slice a row from a matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'r')>
inline void slice_row(Svector<T> &slice, const Tmat &a, Long_I row)
{
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= a.nrows())
		SLS_ERR("out of bound!");
#endif
	Long Nc = a.ncols();
	slice.set(a.ptr() + row * Nc, Nc);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'r')>
inline Svector<T> slice_row(const Tmat &a, Long_I row)
{
	Svector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a col from a matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
inline void slice_col(Svector<T> &slice, const Tmat &a, Long_I col)
{
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= a.ncols())
		SLS_ERR("out of bound!");
#endif
	Long Nr = a.nrows();
	slice.set(a.ptr() + col * Nr, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
inline Svector<T> slice_col(const Tmat &a, Long_I col)
{
	Svector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a Dmat from a mat
// only works for column major for now
template <class Tsmat, class Tmat, SLS_IF(
	is_slice_mat<Tsmat>() && major<Tsmat>() == 'c' &&
	is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
inline void slice_mat(Tsmat &slice, const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Tsmat slice_mat(&a(i, j), Nr, Nc, a.nrows());
}

// slice a3(:,:,i3)
// only supports Cmat<> and Cmat3<> for now
template <class Tmat, class Tmat3, SLS_IF(
	is_same<contain_type<Tmat>, contain_type<Tmat3>>() &&
	is_Scmat<Tmat>() && is_Cmat3d<Tmat3>())>
inline void slice_mat12(Tmat &a, const Tmat3 &a3, Long_I i3)
{
#ifdef SLS_CHECK_BOUNDS
	if (i3 < 0 || i3 >= a.dim3())
		SLS_ERR("out of bound!");
#endif
	a.set(a.ptr(), a.dim1(), a.dim2());
}

} // namespace slisc
