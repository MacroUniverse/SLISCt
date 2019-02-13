// sparse matrix containers
#pragma once
#include "vector.h"

#ifndef NDEBUG
// make sure (i,j) element doesn't exist when using MatCoo<T>::push(s,i,j) 
#define SLS_CHECK_COO_REPEAT
#endif

namespace slisc {

// square diagonal matrix
// mostly a clone a Vector<T>
template <class T>
class Diag : public Vector<T>
{
private:
	typedef Vector<T> Base;
public:
	using Base::operator();
	using Base::operator=;
	Diag() : Base() {}
	Diag(Long_I N) : Base(N) {}
	Diag(Long_I N, const T &s) : Base(N, s) {}
	Diag(Long_I Nr, Long_I Nc) : Base(Nr)
	{
#ifdef SLS_CHECK_SHAPE
		if (Nr != Nc) error("must be a square matrix!");
#endif
	}
	Diag(Long_I Nr, Long_I Nc, const T &s) : Diag(Nr, Nc)
	{ *this = s; }
	Diag(const Vector<T> &v) { *this = v; }
	Long size() const
	{
		error("use nnz() instead!");
	}
	Long nnz() const { return Base::size(); }
	Long nrows() const { return Base::size(); }
	Long ncols() const { return Base::size(); }
	Diag &operator=(const Diag &rhs)
	{ Base::operator=(rhs); return *this; }
	Diag &operator=(const Vector<T> &rhs)
	{ Base::operator=(rhs); return *this; }
	T &operator()(Long_I i, Long_I j)
	{
		if (i == j) return (*this)[i];
		return T();
	}
	const T &operator()(Long_I i, Long_I j) const
	{
		if (i == j) return (*this)[i];
		return T();
	}
};

// convert vector to diagonal matrix
template <class T>
const Diag<T> &diag(const Vector<T> &v)
{ return (Diag<T>&)v; }

// COO sparse matrix
template <class T>
class MatCoo : public Vbase<T>
{
private:
	typedef Vbase<T> Base;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc, m_Nnz;
	Vector<Long> m_row, m_col;
	T m_zero; // TODO: this could be static inline variable for c++17
public:
	using Base::ptr;
	MatCoo(): m_Nr(0), m_Nc(0), m_Nnz(0), m_zero(T()) {}
	MatCoo(Long_I Nr, Long_I Nc) : m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_zero(T()) {}
	MatCoo(Long_I Nr, Long_I Nc, Long_I Nnz):
		Base(Nnz), m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(Nnz), m_col(Nnz), m_zero(T()) {}
	MatCoo(const MatCoo &rhs);		// Copy constructor
	const Long *row_ptr() const { return m_row.ptr(); }
	const Long *col_ptr() const { return m_col.ptr(); }
	MatCoo & operator=(const MatCoo &rhs);
	template <class T1>
	MatCoo & operator=(const MatCoo<T1> &rhs);	// copy assignment (do resize(rhs))
	// inline void operator<<(MatCoo &rhs); // move data and rhs.resize(0, 0); rhs.resize(0)
	T& ref(Long_I i, Long_I j);	// reference to an element
	const T& operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	Long nrows() const { return m_Nr; }
	Long ncols() const { return m_Nc; }
	Long size() const {
		error("use nnz() or capacity() instead!");
	}
	Long nnz() const { return m_Nnz; } // return number of non-zero elements
	Long capacity() const { return Base::size(); }
	T &operator()(Long_I ind); // return element
	const T &operator()(Long_I ind) const;
	Long row(Long_I ind) const; // row index
	Long col(Long_I ind) const; // column index
	void trim(Long_I Nnz); // decrease m_Nnz to Nnz
	void resize(Long_I N)
	{
		error("use reserve instead!");
	}
	void reserve(Long_I N) // reallocate memory, data will be lost
	{ Base::resize(N); m_row.resize(N); m_col.resize(N); m_Nnz = 0; }
	void reshape(Long_I Nr, Long_I Nc) // change matrix shape
	{ m_Nr = Nr; m_Nc = Nc; }
	template <class T1>
	void reserve(const MatCoo<T1> &a);
	template <class T1>
	void reshape(const MatCoo<T1> &a);
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
	return operator=<T>(rhs);
}

template <class T> template <class T1>
inline MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T1> &rhs)
{
	if (this == &rhs) error("self assignment is forbidden!");
	reshape(rhs); reserve(rhs);
	m_row = rhs.m_row;
	m_col = rhs.m_col;
	m_Nnz = rhs.m_Nnz;
	veccpy(ptr(), rhs.ptr(), m_Nnz);
	return *this;
}

template <class T>
inline T& MatCoo<T>::ref(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		error("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n;
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j)
			return m_p[n];
	}
	error("MatCoo::operator()(i,j): element does not exist!");
	return m_p[0];
}

template <class T>
inline const T &MatCoo<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		error("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n;
	for (n = 0; n < m_Nnz; ++n)
		if (row(n) == i && col(n) == j)
			return m_p[n];
	// never return a (const) reference to a temporary
	return m_zero;
}

template <class T>
inline void MatCoo<T>::push(const T &s, Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
		error("MatCoo::push(): index out of bounds!");
#endif
#ifdef SLS_CHECK_COO_REPEAT
	Long n;
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j)
			error("MatCoo::push(s,i,j): element already exists!");
	}
#endif
	if (m_Nnz == m_N) error("MatCoo::add(): out of memory, please reserve!");
	m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
	++m_Nnz;
}

template <class T>
inline void MatCoo<T>::set(const T &s, Long_I i, Long_I j)
{
	Long n;
	// change
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j) {
			m_p[n] = s; m_row[n] = i; m_col[n] = j; return;
		}
	}
	// push
	if (m_Nnz == m_N) error("MatCoo::add(): out of memory, please reserve!");
	m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
	++m_Nnz;
}

template <class T>
inline T & MatCoo<T>::operator()(Long_I ind)
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::operator(): subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
inline const T & MatCoo<T>::operator()(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::operator() const: subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
inline Long MatCoo<T>::row(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::row() subscript out of bounds");
#endif
	return m_row[ind];
}

template <class T>
inline Long MatCoo<T>::col(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind < 0 || ind >= m_Nnz)
		error("MatCoo::col() subscript out of bounds");
#endif
	return m_col[ind];
}

template <class T>
inline void MatCoo<T>::trim(Long_I Nnz)
{
#ifdef SLS_CHECK_SHAPE
	if (Nnz < 0)
		error("MatCoo::trim() negative input!");
#endif
	if (Nnz < m_Nnz) m_Nnz = Nnz;
	else if (Nnz > m_Nnz) error("MatCoo::trim(): Nnz > m_Nnz!");
}

template <class T>
template <class T1>
inline void MatCoo<T>::reshape(const MatCoo<T1> &a)
{
	reshape(a.nrows(), a.ncols());
}

template <class T>
template <class T1>
inline void MatCoo<T>::reserve(const MatCoo<T1> &a)
{
	reserve(a.capacity());
}

// sparse Hermitian / symmetric
// only stores the upper triangle
// nnz() is the actual # of non-zero elem. stored
template <class T>
class MatCooH : public MatCoo<T>
{
private:
	typedef MatCoo<T> Base;
public:
	MatCooH(): Base() {}
	MatCooH(Long_I Nr, Long_I Nc);
	MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz);
	using Base::operator();
	T &ref(Long_I i, Long_I j); // reference to an element
	const T operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	void reshape(Long_I Nr, Long_I Nc)
	{
#ifdef SLS_CHECK_SHAPE
		if (Nr != Nc) error("must be a square matrix!");
#endif
		Base::reshape(Nr, Nc);
	} // change matrix shape
	template <class T1>
	void reshape(const MatCoo<T1> &a);
	template <class T1>
	MatCooH &operator=(const MatCooH<T1> &rhs)
	{ Base(*this) = Base(rhs); }
};

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc) : Base(Nr, Nc)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) error("must be square matrix!");
#endif
}

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz) : Base(Nr, Nc, Nnz)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) error("must be square matrix!");
#endif
}

// cannot return a const reference since conj() might create a temporary
template <class T>
inline const T MatCooH<T>::operator()(Long_I i, Long_I j) const
{
	if (i > j) {
		if constexpr (is_comp<T>())
			return conj(Base::operator()(j, i));
		else
			return Base::operator()(j, i);
	}		
	else
		return Base::operator()(i, j);
}

template <class T>
inline T &MatCooH<T>::ref(Long_I i, Long_I j)
{
	if (i > j)
		error("lower triangle is empty!");
	else
		return Base::ref(i, j);
}

template <class T>
void MatCooH<T>::push(const T &s, Long_I i, Long_I j)
{
	if (i > j)
		Base::push(conj(s), j, i);
	else
		Base::push(s, i, j);
}

template <class T>
void MatCooH<T>::set(const T &s, Long_I i, Long_I j)
{
	if (i > j)
		Base::set(s, j, i);
	else
		Base::set(s, i, j);
}

template <class T> template <class T1>
void MatCooH<T>::reshape(const MatCoo<T1> &a)
{
#ifdef SLS_CHECK_SHAPE
	if (a.nrows() != a.ncols())
		error("a is not square matrix!");
#endif
	reshape(a.nrows());
}

} // namespace slisc
