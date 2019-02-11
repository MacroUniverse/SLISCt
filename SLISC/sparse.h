// sparse matrix classes
#pragma once

#include "arithmetic.h"
#ifndef NDEBUG_
// make sure (i,j) element doesn't exist when using MatCoo<T>::push(s,i,j) 
#define _CHECK_COO_REPEAT_
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
	using Base::size;
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
	Long nrows() const { return size(); }
	Long ncols() const { return size(); }
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

// ptr arithmetics

template <class T, class T1, class T2, SLS_IF(
	is_scalar<T1>() && is_scalar<T2>() &&
	is_same<T, promo_type<T1, T2>>()
)>
void mul_cmat_cmat_diag(T *c, const T1 *a, Long_I Nr, Long_I Nc, const T2 *b)
{
	for (Long i = 0; i < Nc; ++i) {
		times_vvs(c, a, b[i], Nr);
		c += Nr; a += Nr;
	}
}

template <class T, class Tx, class Ty, SLS_IF(
	is_scalar<T>() && is_scalar<Tx>() &&
	is_same<Ty, promo_type<T,Tx>>()
)>
void mul_v_coo_v(Ty *y, const Tx *x, const T *a_ij, const Long *i, const Long *j, Long_I Nr, Long_I N)
{
	vecset(y, Ty(), Nr);
	for (Long k = 0; k < N; ++k)
		y[i[k]] += a_ij[k] * x[j[k]];
}

template <class T, class Tx, class Ty, SLS_IF(
	is_scalar<T>() && is_scalar<Tx>() &&
	is_same<Ty, promo_type<T, Tx>>()
)>
void mul_v_cooh_v(Ty *y, const Tx *x, const T *a_ij, const Long *i, const Long *j, Long_I Nr, Long_I N)
{
	vecset(y, Ty(), Nr);
	for (Long k = 0; k < N; ++k) {
		Long r = i[k], c = j[k];
		if (r == c)
			y[r] += a_ij[k] * x[c];
		else {
			y[r] += a_ij[k] * x[c];
			y[c] += conj(a_ij[k]) * x[r];
		}
	}
}

// arithmetics

template <class T, class Ts, SLS_IF(
	is_MatCoo<T>() && is_scalar<Ts>()
)>
inline void operator*=(T &v, const Ts &s)
{
	times_equals_vs(v.ptr(), s, v.size());
}

// dense matrix +=,-= MatCoo<>

template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_MatCoo<T1>() &&
	is_promo<contain_type<T>, contain_type<T1>>()
)>
inline void operator+=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	for (Long i = 0; i < v1.size(); ++i) {
		v(v1.row(i), v1.col(i)) += v1[i];
	}
}

template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_MatCoo<T1>() &&
	is_promo<contain_type<T>, contain_type<T1>>()
)>
inline void operator-=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1)) error("wrong shape!");
#endif
	for (Long i = 0; i < v1.size(); ++i) {
		v(v1.row(i), v1.col(i)) -= v1[i];
	}
}

// infinite norm
template <class T, SLS_IF(
	type_num<contain_num<T>>() >= 20
)>
inline rm_comp<T> norm_inf(const MatCoo<T> &A)
{
	Vector<rm_comp<T>> abs_sum(A.nrows(), 0);
	for (Long i = 0; i < A.nnz(); ++i) {
		abs_sum(A.row(i)) += abs(A[i]);
	}
	return max(abs_sum);
}

// matrix vector multiplication

template <class Ta, class Tx, class Ty, SLS_IF(
	is_Vector<Ty>() && is_MatCoo<Ta>() && is_Vector<Tx>()
)>
void mul(Ty &y, const Ta &a, const Tx &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != x.size()) error("wrong shape!");
#endif
	y.resize(a.nrows());
	mul_v_coo_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.size());
}

template <class Ta, class Tx, class Ty, SLS_IF(
	is_Vector<Ty>() && is_MatCooH<Ta>() && is_Vector<Tx>()
)>
void mul(Ty &y, const Ta &a, const Tx &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != x.size()) error("wrong shape!");
#endif
	y.resize(a.nrows());
	mul_v_cooh_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.size());
}

// matrix matrix multiplication

// mul(Cmat, Cmat, Diag)
template <class T, class T1, class T2, SLS_IF(
	is_scalar<T>() && is_scalar<T1>() && is_scalar<T2>()
)>
void mul(Cmat<T> &c, const Cmat<T1> &a, const Diag<T2> &b)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != b.nrows()) error("illegal shape!");
#endif
	c.resize(a);
	mul_cmat_cmat_diag(c.ptr(), a.ptr(), a.nrows(), a.ncols(), b.ptr());
}

} // namespace slisc
