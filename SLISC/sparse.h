#pragma once
// sparse matrix classes

#include "slisc.h"
#ifndef NDEBUG_
// make sure (i,j) element doesn't exist when using MatCoo<T>::push(s,i,j) 
#define _CHECK_COO_REPEAT_
#endif

namespace slisc {

template <class T>
class MatCoo : public Vbase<T>
{
private:
	typedef Vbase<T> Base;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc, m_Nnz;
	Vector<Long> m_row, m_col;
public:
	MatCoo(): m_Nr(0), m_Nc(0), m_Nnz(0) {}
	MatCoo(Long_I Nr, Long_I Nc, Long_I N): m_Nr(Nr), m_Nc(Nc), m_Nnz(0), Base(N), m_row(N), m_col(N) {}
	MatCoo(const MatCoo &rhs);		// Copy constructor
	inline MatCoo & operator=(const MatCoo &rhs);	// copy assignment (do resize(rhs))
	// inline void operator<<(MatCoo &rhs); // move data and rhs.resize(0, 0); rhs.resize(0)
	inline T& operator()(Long_I i, Long_I j);	// double indexing (element must exist)
	inline const T& operator()(Long_I i, Long_I j) const;
	inline void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	using Base::size; // return m_N
	Long nrows() const { return m_Nr; }
	Long ncols() const { return m_Nc; }
	Long nnz() const { return m_Nnz; }
	inline T &operator()(Long_I ind); // return element
	inline const T &operator()(Long_I ind) const;
	inline Long & row(Long_I ind); // row index
	inline Long row(Long_I ind) const;
	inline Long & col(Long_I ind); // column index
	inline Long col(Long_I ind) const;
	inline void trim(Long_I Nnz); // decrease m_Nnz to Nnz
	inline void resize(Long_I N) { Base::resize(N); m_row.resize(N); m_col.resize(N); m_Nnz = 0; }
	inline void resize(Long_I Nr, Long_I Nc) // resize (contents preserved)
	{ m_Nr = Nr; m_Nc = Nc; }
	template <class T1>
	inline void resize(const MatCoo<T1> &a);
	~MatCoo() {};
};

template <class T>
MatCoo<T>::MatCoo(const MatCoo<T> &rhs)
{
	error("Copy constructor or move constructor is forbidden, use reference "
		 "argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!");
	resize(rhs);
	m_Nnz = rhs.m_Nnz;
	memcpy(m_p, rhs.m_p, m_Nnz*sizeof(T));
	memcpy(m_row.ptr(), rhs.m_row.ptr(), m_Nnz*sizeof(Long));
	memcpy(m_col.ptr(), rhs.m_col.ptr(), m_Nnz*sizeof(Long));
	return *this;
}

template <class T>
inline T& MatCoo<T>::operator()(Long_I i, Long_I j)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
		error("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n;
	for (n = 0; n < m_Nnz; ++n)
		if (row(n) == i && col(n) == j)
			return m_p[n];
	error("MatCoo::operator()(i,j): element does not exist!");
}

template <class T>
inline const T& MatCoo<T>::operator()(Long_I i, Long_I j) const
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
		error("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n;
	for (n = 0; n < m_Nnz; ++n)
		if (row(n) == i && col(n) == j)
			return m_p[n];
	error("MatCoo::operator()(i,j): element does not exist!");
}

template <class T>
inline void MatCoo<T>::push(const T &s, Long_I i, Long_I j)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
		error("MatCoo::push(): index out of bounds!");
#endif
#ifdef _CHECK_COO_REPEAT_
	Long n;
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j)
			error("MatCoo::push(s,i,j): element already exists!");
	}
#endif
	if (m_Nnz == m_N) error("MatCoo::add(): out of memory, please resize!");
	m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
	++m_Nnz;
}

template <class T>
inline T & MatCoo<T>::operator()(Long_I ind)
{
#ifdef _CHECKBOUNDS_
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::operator(): subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
inline const T & MatCoo<T>::operator()(Long_I ind) const
{
#ifdef _CHECKBOUNDS_
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::operator() const: subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
inline Long & MatCoo<T>::row(Long_I ind)
{
#ifdef _CHECKBOUNDS_
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::row() subscript out of bounds");
#endif
	return m_row(ind);
}

template <class T>
inline Long MatCoo<T>::row(Long_I ind) const
{
#ifdef _CHECKBOUNDS_
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::row() subscript out of bounds");
#endif
	return m_row(ind);
}

template <class T>
inline Long & MatCoo<T>::col(Long_I ind)
{
#ifdef _CHECKBOUNDS_
	if (ind < 0 || ind >= m_Nnz)
		error("MatCoo::col() subscript out of bounds");
#endif
	return m_col(ind);
}

template <class T>
inline Long MatCoo<T>::col(Long_I ind) const
{
#ifdef _CHECKBOUNDS_
	if (ind < 0 || ind >= m_Nnz)
		error("MatCoo::col() subscript out of bounds");
#endif
	return m_col(ind);
}

template <class T>
inline void MatCoo<T>::trim(Long_I Nnz)
{
#ifdef _CHECKBOUNDS_
	if (Nnz < 0)
		error("MatCoo::trim() negative input!");
#endif
	if (Nnz < m_Nnz) m_Nnz = Nnz;
	else if (Nnz > m_Nnz) error("MatCoo::trim(): Nnz > m_Nnz!");
}

template <class T>
template <class T1>
inline void MatCoo<T>::resize(const MatCoo<T1> &a)
{
	resize(a.size());
	resize(a.nrows(), a.ncols()); 
	m_Nnz = 0;
}

typedef MatCoo<Doub> McooDoub, McooDoub_O, McooDoub_IO;
typedef const MatCoo<Doub> McooDoub_I;

// arithmatics

// matrix vector multiplication
void mul(VecDoub_O &v, const McooDoub_I &a, const VecDoub_I &v1)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != v1.size()) error("wrong shape!");
#endif
	Long i;
	v.resize(a.nrows());
	v = 0.;
	for (i = 0; i < a.nnz(); ++i) {
		v(a.row(i)) += a(i) * v1(a.col(i));
	}
}

} // namespace slisc
