#pragma once
#include "vector.h"

namespace slisc {

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
} // namespace slisc
