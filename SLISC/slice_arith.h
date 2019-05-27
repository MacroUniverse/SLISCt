// slicing arithmetics
#pragma once
#include "svector.h"
#include "dvector.h"
#include "scmat.h"
#include "Dcmat.h"

namespace slisc {

// slice a row from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
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
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
inline Svector<T> slice_row(const Tmat &a, Long_I row)
{
	Svector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a row from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
inline void slice_row(Dvector<T> &slice, const Tmat &a, Long_I row)
{
	Long Nr = a.nrows(), Nc = a.ncols();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	slice.set(a.ptr() + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
inline Dvector<T> slice_row(const Tmat &a, Long_I row)
{
	Dvector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a row from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline void slice_row(Dvector<T> &slice, const Tmat &a, Long_I row, Long_I i3 = 0)
{
	Long Nr = a.dim1(), Nc = a.dim2();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	slice.set(a.ptr() + Nr*Nc*i3 + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline Dvector<T> slice_row(const Tmat &a, Long_I row, Long_I i3 = 0)
{
	Dvector<T> slice;
	slice_row(slice, a, row, i3);
	return slice;
}

// slice a col from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
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
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
inline Svector<T> slice_col(const Tmat &a, Long_I col)
{
	Svector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a col from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
inline void slice_col(Dvector<T> &slice, const Tmat &a, Long_I col)
{
	Long Nr = a.nrows(), Nc = a.ncols();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	slice.set(a.ptr() + col, Nc, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
inline Dvector<T> slice_col(const Tmat &a, Long_I col)
{
	Dvector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a col from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline void slice_col(Svector<T> &slice, const Tmat &a, Long_I col, Long_I i3 = 0)
{
	Long Nr = a.dim1(), Nc = a.dim2();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	slice.set(a.ptr() + Nr*Nc*i3 + Nr*col, Nc);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline Svector<T> slice_col(const Tmat &a, Long_I col, Long_I i3 = 0)
{
	Svector<T> slice;
	slice_col(slice, a, col, i3);
	return slice;
}

// slice a Dmat from a mat
// only works for column major for now
template <class Tsmat, class Tmat, SLS_IF(
	is_slice_mat<Tsmat>() && is_cmajor<Tsmat>() &&
	is_dense_mat<Tmat>() && is_cmajor<Tmat>() &&
	is_same_contain_type<Tsmat, Tmat>())>
inline void slice_mat(Tsmat &slice, const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Tsmat slice_mat(&a(i, j), Nr, Nc, a.nrows());
}

// slice a3(:,:,i3)
// only supports Cmat<> and Cmat3<> for now
template <class T, SLS_IF(is_scalar<T>())>
inline void slice_mat12(Scmat<T> &a, const Cmat3d<T> &a3, Long_I i3)
{
#ifdef SLS_CHECK_BOUNDS
	if (i3 < 0 || i3 >= a3.dim3())
		SLS_ERR("out of bound!");
#endif
	a.set_size(a3.dim1(), a3.dim2());
	a.set_ptr(a3.ptr() + a.size() * i3);
}

template <class T, SLS_IF(is_scalar<T>())>
inline Scmat<T> slice_mat12(const Cmat3d<T> &a3, Long_I i3)
{
	Scmat<T> slice;
	slice_mat12(slice, a3, i3);
	return slice;
}

} // namespace slisc
