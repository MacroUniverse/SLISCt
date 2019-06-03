// slicing arithmetics
#pragma once
#include "svector.h"
#include "dvector.h"
#include "scmat.h"
#include "Dcmat.h"

namespace slisc {

// slice a vector from any dense container using single index
template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
void slice_vec(Svector<T> &slice, const Tdens &a, Long_I start, Long_I n)
{
	slice.set(a.ptr() + start, n);
}

template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
Svector<T> slice_vec(const Tdens &a, Long_I start, Long_I n)
{
	Svector<T> slice;
	slice_vec(slice, a, start, n);
	return slice;
}

// slice a Dvector<> from any dense container using single index
template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
void slice_vec(Dvector<T> &slice, const Tdens &a, Long_I start, Long_I n, Long_I step)
{
	slice.set(a.ptr() + start, n, step);
}

template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
Dvector<T> slice_vec(const Tdens &a, Long_I start, Long_I n, Long_I step)
{
	Dvector<T> slice;
	slice_vec(slice, a, start, n, step);
	return slice;
}

// slice a row from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
void slice_row(Svector<T> &slice, const Tmat &a, Long_I row)
{
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= a.n1())
		SLS_ERR("out of bound!");
#endif
	Long Nc = a.n2();
	slice.set(a.ptr() + row * Nc, Nc);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
Svector<T> slice_row(const Tmat &a, Long_I row)
{
	Svector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a row from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice_row(Dvector<T> &slice, const Tmat &a, Long_I row)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	slice.set(a.ptr() + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
Dvector<T> slice_row(const Tmat &a, Long_I row)
{
	Dvector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a row from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
void slice_row(Dvector<T> &slice, const Tmat &a, Long_I row, Long_I k = 0)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	slice.set(a.ptr() + Nr*Nc*k + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
Dvector<T> slice_row(const Tmat &a, Long_I row, Long_I k = 0)
{
	Dvector<T> slice;
	slice_row(slice, a, row, k);
	return slice;
}

// slice a col from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice_col(Svector<T> &slice, const Tmat &a, Long_I col)
{
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= a.n2())
		SLS_ERR("out of bound!");
#endif
	Long Nr = a.n1();
	slice.set(a.ptr() + col * Nr, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
Svector<T> slice_col(const Tmat &a, Long_I col)
{
	Svector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a col from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
void slice_col(Dvector<T> &slice, const Tmat &a, Long_I col)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	slice.set(a.ptr() + col, Nc, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
Dvector<T> slice_col(const Tmat &a, Long_I col)
{
	Dvector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a col from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
void slice_col(Svector<T> &slice, const Tmat &a, Long_I col, Long_I k = 0)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	slice.set(a.ptr() + Nr*Nc*k + Nr*col, Nr);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
Svector<T> slice_col(const Tmat &a, Long_I col, Long_I k = 0)
{
	Svector<T> slice;
	slice_col(slice, a, col, k);
	return slice;
}

// slice a Dmat from a mat
// only works for column major for now
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_cmajor<Tmat>() &&is_dense_mat<Tmat>())>
void slice_mat(Dcmat<T> &slice, const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	slice.set(&a(i, j), Nr, Nc, a.n1());
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_cmajor<Tmat>() &&is_dense_mat<Tmat>())>
Dcmat<T> slice_mat(const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Dcmat<T> slice;
	slice_mat(slice, a, i, Nr, j, Nc);
	return slice;
}

// slice a3(i,j,:)
template <class Tmat3, class T = contain_type<Tmat3>,
	SLS_IF(is_Cmat3d<Tmat3>())>
void slice_dim3(Dvector<T> &slice, const Tmat3 &a,
	Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= a.n1() || j < 0 || j >= a.n2())
		SLS_ERR("index out of bound!");
#endif
	Long N1N2 = a.n1()*a.n2();
	slice.set(a.ptr() + i + a.n1()*j, a.n3(), a.n1()*a.n2());
}

template <class Tmat3, class T = contain_type<Tmat3>,
	SLS_IF(is_Cmat3d<Tmat3>())>
Dvector<T> slice_dim3(const Tmat3 &a,
		Long_I i, Long_I j)
{
	Dvector<T> slice;
	slice_dim3(slice, a, i, j);
	return slice;
}

// slice a3(:,:,k)
// only supports Cmat<> and Cmat3<> for now
template <class T, SLS_IF(is_scalar<T>())>
void slice_mat12(Scmat<T> &a, const Cmat3d<T> &a3, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
	if (k < 0 || k >= a3.n3())
		SLS_ERR("out of bound!");
#endif
	a.set_size(a3.n1(), a3.n2());
	a.set_ptr(a3.ptr() + a.size() * k);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice_mat12(const Cmat3d<T> &a3, Long_I k)
{
	Scmat<T> slice;
	slice_mat12(slice, a3, k);
	return slice;
}

// slice a4(:,:,k,i4)
// only supports Cmat<> and Cmat3<> for now
template <class T, SLS_IF(is_scalar<T>())>
void slice_mat12(Scmat<T> &a, const Cmat4d<T> &a4, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
	if (k < 0 || k >= a4.n3())
		SLS_ERR("out of bound!");
#endif
	a.set_size(a4.n1(), a4.n2());
	a.set_ptr(a4.ptr() + a.size() * (a4.n3()*l  + k));
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice_mat12(const Cmat4d<T> &a4, Long_I k, Long_I l)
{
	Scmat<T> slice;
	slice_mat12(slice, a4, k, l);
	return slice;
}

} // namespace slisc
