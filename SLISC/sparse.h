// sparse matrix classes
#pragma once

#include "arithmetic.h"
#ifndef NDEBUG_
// make sure (i,j) element doesn't exist when using MatCoo<T>::push(s,i,j) 
#define _CHECK_COO_REPEAT_
#endif

namespace slisc {

template <class T> class MatCoo;
typedef MatCoo<Doub> McooDoub;
typedef const McooDoub &McooDoub_I;
typedef McooDoub &McooDoub_O, &McooDoub_IO;

typedef MatCoo<Comp> McooComp;
typedef const McooComp &McooComp_I;
typedef McooComp &McooComp_O, &McooComp_IO;

template <class T> class MatCooH;
typedef MatCooH<Doub> McoohDoub;
typedef const McoohDoub &McoohDoub_I;
typedef McoohDoub &McoohDoub_O, &McoohDoub_IO;

typedef MatCooH<Comp> McoohComp;
typedef const McoohComp &McoohComp_I;
typedef McoohComp &McoohComp_O, &McoohComp_IO;

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
	Diag(Long_I N, const T &s) : Base(N, s) {}
	Diag(const Vector<T> &v) { *this = v; }
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

// Cmat times Diag
template <class T, class T1, class T2>
void mul(Cmat<T> &v, const Cmat<T1> &v1, const Diag<T2> &v2)
{
	Long Nr = v1.nrows(), Nc = v1.ncols();
#ifdef _CHECKBOUNDS_
	if (Nc != v2.size()) error("illegal shape!");
#endif
	v.resize(Nr, v2.size());
	T * p = v.ptr();
	const T1 *p1 = v1.ptr();
	for (Long i = 0; i < Nc; ++i) {
		times_vvs(p, p1, v2[i], Nr);
		p += Nr; p1 += Nr;
	}
}

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
	enum { NDIMS = 2 };
	MatCoo(): m_Nr(0), m_Nc(0), m_Nnz(0), m_zero(T()) {}
	MatCoo(Long_I Nr, Long_I Nc) : m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_zero(T()) {}
	MatCoo(Long_I Nr, Long_I Nc, Long_I Nnz):
		Base(Nnz), m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(Nnz), m_col(Nnz), m_zero(T()) {}
	MatCoo(const MatCoo &rhs);		// Copy constructor
	MatCoo & operator=(const MatCoo &rhs);
	template <class T1>
	MatCoo & operator=(const MatCoo<T1> &rhs);	// copy assignment (do resize(rhs))
	// inline void operator<<(MatCoo &rhs); // move data and rhs.resize(0, 0); rhs.resize(0)
	T& ref(Long_I i, Long_I j);	// reference to an element
	const T& operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	using Base::size; // return m_N
	Long nrows() const { return m_Nr; }
	Long ncols() const { return m_Nc; }
	Long nnz() const { return m_Nnz; }
	T &operator()(Long_I ind); // return element
	const T &operator()(Long_I ind) const;
	Long row(Long_I ind) const; // row index
	Long col(Long_I ind) const; // column index
	void trim(Long_I Nnz); // decrease m_Nnz to Nnz
	void resize(Long_I N) // reallocate memory
	{ Base::resize(N); m_row.resize(N); m_col.resize(N); m_Nnz = 0; }
	void reshape(Long_I Nr, Long_I Nc) // change matrix size
	{ m_Nr = Nr; m_Nc = Nc; }
	template <class T1>
	void resize(const MatCoo<T1> &a);
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
	reshape(rhs); resize(rhs);
	m_row = rhs.m_row;
	m_col = rhs.m_col;
	m_Nnz = rhs.m_Nnz;
	for (Long i = 0; i < m_Nnz; ++i) {
		(*this)(i) = rhs(i);
	}
	return *this;
}

template <class T>
inline T& MatCoo<T>::ref(Long_I i, Long_I j)
{
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
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
#ifdef _CHECKBOUNDS_
	if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
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
inline Long MatCoo<T>::row(Long_I ind) const
{
#ifdef _CHECKBOUNDS_
	if (ind<0 || ind>=m_Nnz)
		error("MatCoo::row() subscript out of bounds");
#endif
	return m_row[ind];
}

template <class T>
inline Long MatCoo<T>::col(Long_I ind) const
{
#ifdef _CHECKBOUNDS_
	if (ind < 0 || ind >= m_Nnz)
		error("MatCoo::col() subscript out of bounds");
#endif
	return m_col[ind];
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
inline void MatCoo<T>::reshape(const MatCoo<T1> &a)
{
	reshape(a.nrows(), a.ncols());
}

template <class T>
template <class T1>
inline void MatCoo<T>::resize(const MatCoo<T1> &a)
{
	resize(a.size());
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
	enum { NDIMS = 2 };
	MatCooH(): Base() {}
	MatCooH(Long_I Nr, Long_I Nc);
	MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz);
	using Base::operator();
	T &ref(Long_I i, Long_I j); // reference to an element
	const T operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
	void push(const T &s, Long_I i, Long_I j); // add one nonzero element
	void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
	void reshape(Long_I N) { Base::reshape(N, N); } // change matrix size
	template <class T1>
	void reshape(const MatCoo<T1> &a);
	template <class T1>
	MatCooH &operator=(const MatCooH<T1> &rhs)
	{ Base(*this) = Base(rhs); }
};

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc) : Base(Nr, Nc)
{
#ifdef _CHECKBOUNDS_
	if (Nr != Nc) error("must be square matrix!");
#endif
}

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz) : Base(Nr, Nc, Nnz)
{
#ifdef _CHECKBOUNDS_
	if (Nr != Nc) error("must be square matrix!");
#endif
}

// cannot return a const reference since conj() might create a temporary
template <class T>
inline const T MatCooH<T>::operator()(Long_I i, Long_I j) const
{
	if (i > j)
		// conj has no overhead for non-complex argument
		return conj(Base::operator()(j, i));
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
#ifdef _CHECKBOUNDS_
	if (a.nrows() != a.ncols())
		error("a is not square matrix!");
#endif
	reshape(a.nrows());
}

// arithmetics

// matrix vector multiplication
// internal only: no bound checking!
template <class T>
void coo_mul(T *y, const MatCoo<T> &a, const T *x)
{
	Long i;
	vecset(y, T(), a.nrows());
	for (i = 0; i < a.nnz(); ++i) {
		y[a.row(i)] += a(i) * x[a.col(i)];
	}
}

// matrix vector multiplication
template <class T, class T1, class T2>
void mul(Vector<T> &y, const MatCoo<T1> &a, const Vector<T2> &x)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != x.size()) error("wrong shape!");
#endif
	y.resize(a.nrows());
	coo_mul(y.ptr(), a, x.ptr());
}

// internal only: no bound checking!
template <class T, class T1, class T2>
void cooh_mul(T *y, const MatCooH<T1> &a, const T2 *x)
{
	Long i;
	vecset(y, T(), a.nrows());
	for (i = 0; i < a.nnz(); ++i) {
		Long r = a.row(i), c = a.col(i);
		if (r == c)
			y[r] += a(i) * x[c];
		else {
			y[r] += a(i) * x[c];
			y[c] += conj(a(i)) * x[r];
		}
	}
}

template <class T, class T1, class T2>
void mul(Vector<T> &y, const MatCooH<T1> &a, const Vector<T2> &x)
{
#ifdef _CHECKBOUNDS_
	if (a.ncols() != x.size()) error("wrong shape!");
#endif
	y.resize(a.nrows());
	cooh_mul(y.ptr(), a, x.ptr());
}

template <class T, class T1>
inline void operator*=(MatCoo<T> &v, const T1 &s)
{ times_equals1(v.ptr(), s, v.size()); }

// dense matrix - sparse matrix
template <class T, class T1>
inline void operator-=(Matrix<T> &v, const MatCoo<T1> &v1)
{
#ifdef _CHECKBOUNDS_
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	for (Long i = 0; i < v1.size(); ++i) {
		v(v1.row(i), v1.col(i)) -= v1(i);
	}
}

inline Doub norm_inf(McooComp_I A)
{
	VecDoub abs_sum(A.nrows(), 0.);
	for (Int i = 1; i < A.nnz(); ++i) {
		abs_sum(A.row(i)) += std::abs(A(i));
	}
	return max(abs_sum);
}

} // namespace slisc
