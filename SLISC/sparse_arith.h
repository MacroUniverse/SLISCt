// arithmetics for sparse containers

#pragma once
#include "sparse.h"
#include "ptr_arith.h"

namespace slisc {

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

template <class T, class T1, class T2, SLS_IF(
	is_scalar<T1>() && is_scalar<T2>() &&
	is_same<T, promo_type<T1, T2>>()
)>
void mul_cmat_diag_cmat(T *c, const T1 *a, const T2 *b, Long_I Nr, Long_I Nc)
{
	for (Long i = 0; i < Nc; ++i) {
		times_vvv(c, b, a, Nr);
		c += Nr; b += Nr;
	}
}

template <class T, class Tx, class Ty, SLS_IF(
	is_scalar<T>() && is_scalar<Tx>() &&
	is_same<Ty, promo_type<T,Tx>>()
)>
void mul_v_coo_v(Ty *y, const Tx *x, const T *a_ij, const Long *i, const Long *j, Long_I Nr, Long_I Nnz)
{
	vecset(y, Ty(), Nr);
	for (Long k = 0; k < Nnz; ++k)
		y[i[k]] += a_ij[k] * x[j[k]];
}

template <class T, class Tx, class Ty, SLS_IF(
	is_scalar<T>() && is_scalar<Tx>() &&
	is_same<Ty, promo_type<T, Tx>>()
)>
void mul_v_cooh_v(Ty *y, const Tx *x, const T *a_ij, const Long *i, const Long *j, Long_I Nr, Long_I Nnz)
{
	vecset(y, Ty(), Nr);
	for (Long k = 0; k < Nnz; ++k) {
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
	times_equals_vs(v.ptr(), s, v.nnz());
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
	type_num<T>() >= 20
)>
inline rm_comp<T> norm_inf(const MatCoo<T> &A)
{
	Vector<rm_comp<T>> abs_sum(A.nrows(), 0.);
	for (Long i = 0; i < A.nnz(); ++i) {
		abs_sum(A.row(i)) += abs(A[i]);
	}
	return max(abs_sum);
}

template <class T, SLS_IF(
	type_num<T>() >= 20
)>
inline rm_comp<T> norm_inf(const MatCooH<T> &A)
{
	Vector<rm_comp<T>> abs_sum(A.nrows(), 0.);
	for (Long i = 0; i < A.nnz(); ++i) {
		Long r = A.row(i), c = A.col(i);
		auto val = abs(A[i]);
		abs_sum(r) += val;
		if (r != c)
			abs_sum(c) += val;
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
	mul_v_coo_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.nnz());
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
	mul_v_cooh_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.nnz());
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

// mul(Cmat, Diag, Cmat)
template <class T, class T1, class T2, SLS_IF(
	is_scalar<T>() && is_scalar<T1>() && is_scalar<T2>()
)>
void mul(Cmat<T> &c, const Diag<T2> &a, const Cmat<T1> &b)
{
#ifdef SLS_CHECK_SHAPE
	if (a.ncols() != b.nrows()) error("illegal shape!");
#endif
	c.resize(b);
	mul_cmat_diag_cmat(c.ptr(), a.ptr(), b.ptr(), b.nrows(), b.ncols());
}

} // namespace slisc
