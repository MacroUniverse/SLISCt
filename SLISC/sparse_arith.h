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

// a(blk_size, blk_size, Nblk) is column major
// overlapped element already divided by 2
template <class T, class Tx, class Ty, SLS_IF(
	is_scalar<T>() && is_scalar<Tx>() &&
	is_same<Ty, promo_type<T, Tx>>()
)>
void mul_v_obd_v(Ty *y, const Tx *x, const T *a, Long_I blk_size, Long_I Nblk)
{
	vecset(y, Ty(0), (blk_size - 1) * Nblk + 1);
	Long step = blk_size - 1;
	for (Long blk = 0; blk < Nblk; ++blk) {
		for (Long j = 0; j < blk_size; ++j) {
			Tx s = x[j];
			for (Long i = 0; i < blk_size; ++i) {
				y[i] += (*a) * s;
				++a;
			}
		}
		x += step; y += step;
	}
}

template <class Tx, class Ty, class Ta, SLS_IF(
	is_dense_vec<Tx>() && is_dense_vec<Ty>() && is_CmatObd<Ta>())>
void mul(Ty &y, const Ta &a, const Tx &x)
{
	mul_v_obd_v(y.ptr(), x.ptr(), a.ptr(), a.n0(), a.nblk());
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
	if (!shape_cmp(v, v1)) SLS_ERR("wrong shape!");
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
	if (!shape_cmp(v, v1)) SLS_ERR("wrong shape!");
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
	Vector<rm_comp<T>> abs_sum(A.n1(), 0.);
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
	Vector<rm_comp<T>> abs_sum(A.n1(), 0.);
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
	is_Vector<Ty>() && is_MatCoo<Ta>() && is_dense_vec<Tx>()
)>
void mul(Ty &y, const Ta &a, const Tx &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != x.size() || a.n1() != y.size())
		SLS_ERR("wrong shape!");
#endif
	mul_v_coo_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

template <class Ta, class Tx, class Ty, SLS_IF(
	is_dense_vec<Ty>() && is_MatCooH<Ta>() && is_dense_vec<Tx>()
)>
void mul(Ty &y, const Ta &a, const Tx &x)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != x.size() || a.rows() != y.size()) SLS_ERR("wrong shape!");
#endif
	mul_v_cooh_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

// matrix matrix multiplication

// mul(Cmat, Cmat, Diag)
template <class T, class T1, class T2, SLS_IF(
	is_scalar<T>() && is_scalar<T1>() && is_scalar<T2>()
)>
void mul(Cmat<T> &c, const Cmat<T1> &a, const Diag<T2> &b)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != b.n1()) SLS_ERR("illegal shape!");
	if (c.n1() != a.n1() || c.n2() != b.n2())
		SLS_ERR("illegal shape!");
#endif
	mul_cmat_cmat_diag(c.ptr(), a.ptr(), a.n1(), a.n2(), b.ptr());
}

// mul(Cmat, Diag, Cmat)
template <class T, class T1, class T2, SLS_IF(
	is_scalar<T>() && is_scalar<T1>() && is_scalar<T2>()
)>
void mul(Cmat<T> &c, const Diag<T2> &a, const Cmat<T1> &b)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n2() != b.n1()) SLS_ERR("illegal shape!");
#endif
	c.resize(b);
	mul_cmat_diag_cmat(c.ptr(), a.ptr(), b.ptr(), b.n1(), b.n2());
}

} // namespace slisc
