// slicing arithmetics
#pragma once
#include "svector.h"
#include "dvector.h"
#include "scmat.h"
#include "Dcmat.h"

namespace slisc {

// slice a row from a matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'r')>
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
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'r')>
inline Svector<T> slice_row(const Tmat &a, Long_I row)
{
	Svector<T> slice;
	slice_row(slice, a, row);
	return slice;
}

// slice a col from a matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
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
	SLS_IF(is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
inline Svector<T> slice_col(const Tmat &a, Long_I col)
{
	Svector<T> slice;
	slice_col(slice, a, col);
	return slice;
}

// slice a Dmat from a mat
// only works for column major for now
template <class Tsmat, class Tmat, SLS_IF(
	is_slice_mat<Tsmat>() && major<Tsmat>() == 'c' &&
	is_dense_mat<Tmat>() && major<Tmat>() == 'c')>
inline void slice_mat(Tsmat &slice, const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Tsmat slice_mat(&a(i, j), Nr, Nc, a.nrows());
}

// slice a3(:,:,i3)
// only supports Cmat<> and Cmat3<> for now
template <class Tmat, class Tmat3, SLS_IF(
	is_same<contain_type<Tmat>, contain_type<Tmat3>>() &&
	is_Scmat<Tmat>() && is_Cmat3d<Tmat3>())>
inline void slice_mat12(Tmat &a, const Tmat3 &a3, Long_I i3)
{
#ifdef SLS_CHECK_BOUNDS
	if (i3 < 0 || i3 >= a.dim3())
		SLS_ERR("out of bound!");
#endif
	a.set(a.ptr(), a.dim1(), a.dim2());
}
} // namespace slisc
