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
	Diag() : Base() {} // default constructor: uninitialized
public:
	using Base::operator();
	using Base::operator=;
	Diag(Long_I N);
	Diag(Long_I N, const T &s);
	Diag(Long_I Nr, Long_I Nc);
	Diag(Long_I Nr, Long_I Nc, const T &s);
	Diag(const Vector<T> &v);
	Long size() const;
	Long nnz() const;
	Long n1() const;
	Long ncols() const;
	Diag &operator=(const Diag &rhs);
	Diag &operator=(const Vector<T> &rhs);
	T &operator()(Long_I i, Long_I j);
	const T &operator()(Long_I i, Long_I j) const;
};

template <class T>
Diag<T>::Diag(Long_I N) : Base(N) {}

template <class T>
Diag<T>::Diag(Long_I N, const T &s) : Base(N, s) {}

template <class T>
Diag<T>::Diag(Long_I Nr, Long_I Nc) : Base(Nr)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) SLS_ERR("must be a square matrix!");
#endif
}

template <class T>
Diag<T>::Diag(Long_I Nr, Long_I Nc, const T &s) : Diag(Nr, Nc)
{
	*this = s;
}

template <class T>
Diag<T>::Diag(const Vector<T> &v) : Base(v.size())
{
	*this = v;
}

template <class T>
Long Diag<T>::size() const
{
	SLS_ERR("use nnz() instead!");
}

template <class T>
Long Diag<T>::nnz() const
{
	return Base::size();
}

template <class T>
Long Diag<T>::n1() const
{
	return Base::size();
}

template <class T>
Long Diag<T>::ncols() const
{
	return Base::size();
}

template <class T>
Diag<T> &Diag<T>::operator=(const Diag &rhs)
{
	Base::operator=(rhs); return *this;
}

template <class T>
Diag<T> &Diag<T>::operator=(const Vector<T> &rhs)
{
	Base::operator=(rhs); return *this;
}

template <class T>
T &Diag<T>::operator()(Long_I i, Long_I j)
{
	if (i == j) return (*this)[i];
	return T();
}

template <class T>
const T &Diag<T>::operator()(Long_I i, Long_I j) const
{
	if (i == j) return (*this)[i];
	return T();
}

// convert vector to diagonal matrix
template <class T>
const Diag<T> &diag(const Vector<T> &v)
{
	return (Diag<T>&)v;
}

// COO sparse matrix
template <class T>
class MatCoo : public Vbase<T>
{
private:
	typedef Vbase<T> Base;
	using Base::m_p;
	using Base::m_N;
	Long m_Nr, m_Nc, m_Nnz;
	VecLong m_row, m_col;
	static inline T m_zero = (T)0; // TODO: this could be static inline variable for c++17
	MatCoo() {} // default constructor: uninitialized
public:
	using Base::ptr;
	MatCoo(Long_I Nr, Long_I Nc);
	MatCoo(Long_I Nr, Long_I Nc, Long_I Nnz);
	MatCoo(const MatCoo &rhs);		// Copy constructor
	const Long *row_ptr() const;
	const Long *col_ptr() const;
	MatCoo & operator=(const MatCoo &rhs);
	template <class T1>
	MatCoo & operator=(const MatCoo<T1> &rhs);	// copy assignment (do resize(rhs))
	// inline void operator<<(MatCoo &rhs); // move data and rhs.resize(0, 0); rhs.resize(0)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	Long n1() const;
	Long ncols() const;
	Long size() const; // forbidden
	Long nnz() const; // return number of non-zero elements
	Long capacity() const;
	// get single index using double index, return -1 if not found
	Long find(Long_I i, Long_I j) const;
	// reference to an element (element must exist)
	T& ref(Long_I i, Long_I j);
	// double indexing (element need not exist)
	const T& operator()(Long_I i, Long_I j) const;
	T &operator()(Long_I ind); // return element
	const T &operator()(Long_I ind) const;
	Long row(Long_I ind) const; // row index
	Long col(Long_I ind) const; // column index
	void trim(Long_I Nnz); // decrease m_Nnz to Nnz
	void resize(Long_I N); // forbidden
	void reserve(Long_I N); // reallocate memory, data will be lost
	void reshape(Long_I Nr, Long_I Nc); // change matrix shape
	template <class T1>
	void reserve(const MatCoo<T1> &a);
	template <class T1>
	void reshape(const MatCoo<T1> &a);
};

template <class T>
MatCoo<T>::MatCoo(Long_I Nr, Long_I Nc)
	: m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(0), m_col(0)
{
	m_N = 0;
}

template <class T>
MatCoo<T>::MatCoo(Long_I Nr, Long_I Nc, Long_I Nnz) :
	Base(Nnz), m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(Nnz), m_col(Nnz) {}

template <class T>
MatCoo<T>::MatCoo(const MatCoo<T> &rhs)
{
	SLS_ERR("Copy constructor or move constructor is forbidden, use reference "
		 "argument for function input or output, and use \"=\" to copy!");
}

template <class T>
const Long *MatCoo<T>::row_ptr() const
{
	return m_row.ptr();
}

template <class T>
const Long *MatCoo<T>::col_ptr() const
{
	return m_col.ptr();
}

template <class T>
MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T> &rhs)
{
	return operator=<T>(rhs);
}

template <class T> template <class T1>
MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T1> &rhs)
{
	if ((void*)this == (void*)&rhs) SLS_ERR("self assignment is forbidden!");
	reshape(rhs); reserve(rhs);
	m_row = rhs.n1();
	m_col = rhs.ncols();
	m_Nnz = rhs.nnz();
	veccpy(ptr(), rhs.ptr(), m_Nnz);
	veccpy(m_row.ptr(), rhs.row_ptr(), m_Nnz);
	veccpy(m_col.ptr(), rhs.col_ptr(), m_Nnz);
	return *this;
}

template <class T>
Long MatCoo<T>::find(Long_I i, Long_I j) const
{
	for (Long n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j)
			return n;
	}
	return -1;
}

template <class T>
T& MatCoo<T>::ref(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n = find(i, j);
	if (n < 0)
		SLS_ERR("MatCoo::operator()(i,j): element does not exist!");
	return m_p[n];
}

template <class T>
const T &MatCoo<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
		SLS_ERR("MatCoo::operator()(i,j): index out of bounds!");
#endif
	Long n = find(i, j);
	if (n < 0)
		return m_zero; // never return a (const) reference to a temporary
	return m_p[n];
}

template <class T>
void MatCoo<T>::push(const T &s, Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
		SLS_ERR("MatCoo::push(): index out of bounds!");
#endif
#ifdef SLS_CHECK_COO_REPEAT
	Long n;
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j)
			SLS_ERR("MatCoo::push(s,i,j): element already exists!");
	}
#endif
	if (m_Nnz == m_N) SLS_ERR("MatCoo::add(): out of memory, please reserve!");
	m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
	++m_Nnz;
}

template <class T>
void MatCoo<T>::set(const T &s, Long_I i, Long_I j)
{
	Long n;
	// change
	for (n = 0; n < m_Nnz; ++n) {
		if (row(n) == i && col(n) == j) {
			m_p[n] = s; m_row[n] = i; m_col[n] = j; return;
		}
	}
	// push
	if (m_Nnz == m_N) SLS_ERR("MatCoo::add(): out of memory, please reserve!");
	m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
	++m_Nnz;
}

template <class T>
Long MatCoo<T>::n1() const
{
	return m_Nr;
}

template <class T>
Long MatCoo<T>::ncols() const
{
	return m_Nc;
}

template <class T>
Long MatCoo<T>::size() const
{
	SLS_ERR("use nnz() or capacity() instead!");
}

template <class T>
Long MatCoo<T>::nnz() const
{
	return m_Nnz;
}

template <class T>
Long MatCoo<T>::capacity() const
{
	return Base::size();
}

template <class T>
T & MatCoo<T>::operator()(Long_I ind)
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		SLS_ERR("MatCoo::operator(): subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
const T & MatCoo<T>::operator()(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		SLS_ERR("MatCoo::operator() const: subscript out of bounds!");
#endif
	return m_p[ind];
}

template <class T>
Long MatCoo<T>::row(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind<0 || ind>=m_Nnz)
		SLS_ERR("MatCoo::row() subscript out of bounds");
#endif
	return m_row[ind];
}

template <class T>
Long MatCoo<T>::col(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
	if (ind < 0 || ind >= m_Nnz)
		SLS_ERR("MatCoo::col() subscript out of bounds");
#endif
	return m_col[ind];
}

template <class T>
void MatCoo<T>::trim(Long_I Nnz)
{
#ifdef SLS_CHECK_SHAPE
	if (Nnz < 0)
		SLS_ERR("MatCoo::trim() negative input!");
#endif
	if (Nnz < m_Nnz) m_Nnz = Nnz;
	else if (Nnz > m_Nnz) SLS_ERR("MatCoo::trim(): Nnz > m_Nnz!");
}

template <class T>
void MatCoo<T>::resize(Long_I N)
{
	SLS_ERR("use reserve instead!");
}

template <class T>
void MatCoo<T>::reserve(Long_I N)
{
	Base::resize(N); m_row.resize(N); m_col.resize(N); m_Nnz = 0;
}

template <class T>
void MatCoo<T>::reshape(Long_I Nr, Long_I Nc)
{
	m_Nr = Nr; m_Nc = Nc;
}

template <class T>
template <class T1>
void MatCoo<T>::reshape(const MatCoo<T1> &a)
{
	reshape(a.n1(), a.ncols());
}

template <class T>
template <class T1>
void MatCoo<T>::reserve(const MatCoo<T1> &a)
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
	MatCooH() : Base() {} // default constructor : uninitialized
public:
	MatCooH(Long_I Nr, Long_I Nc);
	MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz);
	using Base::operator();
	T &ref(Long_I i, Long_I j); // reference to an element
	const T operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	void reshape(Long_I Nr, Long_I Nc); // change matrix shape
	template <class T1>
	void reshape(const MatCoo<T1> &a);
	template <class T1>
	MatCooH &operator=(const MatCooH<T1> &rhs);
};

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc) : Base(Nr, Nc)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) SLS_ERR("must be square matrix!");
#endif
}

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz) : Base(Nr, Nc, Nnz)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) SLS_ERR("must be square matrix!");
#endif
}

// cannot return a const reference since conj() might create a temporary
template <class T>
const T MatCooH<T>::operator()(Long_I i, Long_I j) const
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
T &MatCooH<T>::ref(Long_I i, Long_I j)
{
	if (i > j)
		SLS_ERR("lower triangle is empty!");
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

template <class T>
void MatCooH<T>::reshape(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
	if (Nr != Nc) SLS_ERR("must be a square matrix!");
#endif
	Base::reshape(Nr, Nc);
}

template <class T> template <class T1>
void MatCooH<T>::reshape(const MatCoo<T1> &a)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n1() != a.ncols())
		SLS_ERR("a is not square matrix!");
#endif
	reshape(a.n1());
}

template <class T> template <class T1>
MatCooH<T> &MatCooH<T>::operator=(const MatCooH<T1> &rhs)
{
	Base(*this) = Base(rhs);
}

} // namespace slisc
