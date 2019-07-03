// slicing arithmetics
#pragma once
#include "svector.h"
#include "dvector.h"
#include "scmat.h"
#include "cmat4d.h"
#include "dcmat.h"
#include "jcmat.h"
#include "jcmat3d.h"
#include "jcmat4d.h"

namespace slisc {

// slice a vector from any dense container using single index
template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
void slice_vec(Svector<T> &sli, const Tdens &a, Long_I start, Long_I n)
{
	sli.set(a.ptr() + start, n);
}

template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
Svector<T> slice_vec(const Tdens &a, Long_I start, Long_I n)
{
	Svector<T> sli;
	slice_vec(sli, a, start, n);
	return sli;
}

// slice a Dvector<> from any dense container using single index
template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
void slice_vec(Dvector<T> &sli, const Tdens &a, Long_I start, Long_I n, Long_I step)
{
	sli.set(a.ptr() + start, n, step);
}

template <class Tdens, class T = contain_type<Tdens>,
	SLS_IF(is_dense<Tdens>())>
Dvector<T> slice_vec(const Tdens &a, Long_I start, Long_I n, Long_I step)
{
	Dvector<T> sli;
	slice_vec(sli, a, start, n, step);
	return sli;
}

// slice a row from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
void slice2(Svector<T> &sli, const Tmat &a, Long_I row)
{
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= a.n1())
		SLS_ERR("out of bound!");
#endif
	Long Nc = a.n2();
	sli.set(a.ptr() + row * Nc, Nc);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
Svector<T> slice2(const Tmat &a, Long_I row)
{
	Svector<T> sli;
	slice2(sli, a, row);
	return sli;
}

// slice a row from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice2(Dvector<T> &sli, const Tmat &a, Long_I row)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	sli.set(a.ptr() + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
Dvector<T> slice2(const Tmat &a, Long_I row)
{
	Dvector<T> sli;
	slice2(sli, a, row);
	return sli;
}

// slice a row from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
void slice2(Dvector<T> &sli, const Tmat &a, Long_I row, Long_I k = 0)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (row < 0 || row >= Nr)
		SLS_ERR("out of bound!");
#endif
	
	sli.set(a.ptr() + Nr*Nc*k + row, Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
Dvector<T> slice2(const Tmat &a, Long_I row, Long_I k = 0)
{
	Dvector<T> sli;
	slice2(sli, a, row, k);
	return sli;
}

// slice a col from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice1(Svector<T> &sli, const Tmat &a, Long_I col)
{
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= a.n2())
		SLS_ERR("out of bound!");
#endif
	Long Nr = a.n1();
	sli.set(a.ptr() + col * Nr, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
Svector<T> slice1(const Tmat &a, Long_I col)
{
	Svector<T> sli;
	slice1(sli, a, col);
	return sli;
}

// slice a col from a row-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
void slice1(Dvector<T> &sli, const Tmat &a, Long_I col)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	sli.set(a.ptr() + col, Nc, Nr);
}

// note that `slice1 = slice2` will copy data instead of copying slice object
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_dense_mat<Tmat>() && is_rmajor<Tmat>())>
Dvector<T> slice1(const Tmat &a, Long_I col)
{
	Dvector<T> sli;
	slice1(sli, a, col);
	return sli;
}

// slice a col from a col-major 3d array
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
void slice1(Svector<T> &sli, const Tmat &a, Long_I col, Long_I k = 0)
{
	Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
	if (col < 0 || col >= Nc)
		SLS_ERR("out of bound!");
#endif
	sli.set(a.ptr() + Nr*Nc*k + Nr*col, Nr);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
	is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
Svector<T> slice1(const Tmat &a, Long_I col, Long_I k)
{
	Svector<T> sli;
	slice1(sli, a, col, k);
	return sli;
}

// slice a Dmat from a dense mat
// only works for column major for now
template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_cmajor<Tmat>() &&is_dense_mat<Tmat>())>
void slice(Dcmat<T> &sli, const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	sli.set(&a(i, j), Nr, Nc, a.n1());
}

template <class Tmat, class T = contain_type<Tmat>,
	SLS_IF(is_cmajor<Tmat>() &&is_dense_mat<Tmat>())>
Dcmat<T> slice(const Tmat &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Dcmat<T> sli;
	slice(sli, a, i, Nr, j, Nc);
	return sli;
}

// slice a Dcmat from a Dcmat
// only works for column major for now
template <class T, SLS_IF(is_scalar<T>())>
void slice(Dcmat<T> &sli, const Dcmat<T> &a,
		Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	sli.set(&a(i, j), Nr, Nc, a.lda());
}

template <class T, SLS_IF(is_scalar<T>())>
Dcmat<T> slice(const Dcmat<T> &a,
		Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Dcmat<T> sli;
	slice(sli, a, i, Nr, j, Nc);
	return sli;
}

// slice a Jcmat from a Jcmat
template <class T, SLS_IF(is_scalar<T>())>
void slice(Jcmat<T> &sli, const Jcmat<T> &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	sli.set(&a(i, j), Nr, Nc, a.step1(), a.step2());
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat<T> slice(const Jcmat<T> &a,
	Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
	Jcmat<T> sli;
	slice(sli, a, i, Nr, j, Nc);
	return sli;
}

// slice a3(i,j,:)
template <class Tmat3, class T = contain_type<Tmat3>,
	SLS_IF(is_Cmat3d<Tmat3>())>
void slice3(Dvector<T> &sli, const Tmat3 &a,
	Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= a.n1() || j < 0 || j >= a.n2())
		SLS_ERR("index out of bound!");
#endif
	Long N1N2 = a.n1()*a.n2();
	sli.set(a.ptr() + i + a.n1()*j, a.n3(), a.n1()*a.n2());
}

template <class Tmat3, class T = contain_type<Tmat3>,
	SLS_IF(is_Cmat3d<Tmat3>())>
Dvector<T> slice3(const Tmat3 &a,
		Long_I i, Long_I j)
{
	Dvector<T> sli;
	slice3(sli, a, i, j);
	return sli;
}

// slice a3(:,:,k)
// only supports Cmat<> and Cmat3<> for now
template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat<T> &a, const Cmat3d<T> &a3, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
	if (k < 0 || k >= a3.n3())
		SLS_ERR("out of bound!");
#endif
	a.set_size(a3.n1(), a3.n2());
	a.set_ptr(a3.ptr() + a.size() * k);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice12(const Cmat3d<T> &a3, Long_I k)
{
	Scmat<T> sli;
	slice12(sli, a3, k);
	return sli;
}

// slice a3(i,:,:)
template <class T, SLS_IF(is_scalar<T>())>
void slice23(Jcmat<T> &sli, const Cmat3d<T> &a3, Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= a3.n1())
		SLS_ERR("out of bound!");
#endif
	sli.set(a3.ptr() + i, a3.n2(), a3.n3(), a3.n1(), a3.n1()*a3.n2());
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat<T> slice23(const Cmat3d<T> &a3, Long_I i)
{
	Jcmat<T> sli;
	slice23(sli, a3, i);
	return sli;
}

// slice a4(i,:,:,l)
template <class T, SLS_IF(is_scalar<T>())>
void slice23(Jcmat<T> &sli, const Cmat4d<T> &a4, Long_I i, Long_I l)
{
	Long N1 = a4.n1(), N2 = a4.n2(), N3 = a4.n3(), N4 = a4.n4();
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= N1 || l < 0 || l > N4)
		SLS_ERR("out of bound!");
#endif
	Long N1N2 = N1 * N2;
	sli.set(a4.ptr() + i + N1N2*N3*l, N2, N3, N1, N1N2);
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat<T> slice23(const Cmat4d<T> &a4, Long_I i, Long_I l)
{
	Jcmat<T> sli;
	slice23(sli, a4, i, l);
	return sli;
}

// slice a4(:,:,k,l)
// only supports Cmat<> and Cmat3<> for now
template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat<T> &a, const Cmat4d<T> &a4, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
	if (k < 0 || k >= a4.n3())
		SLS_ERR("out of bound!");
#endif
	a.set_size(a4.n1(), a4.n2());
	a.set_ptr(a4.ptr() + a.size() * (a4.n3()*l  + k));
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice12(const Cmat4d<T> &a4, Long_I k, Long_I l)
{
	Scmat<T> sli;
	slice12(sli, a4, k, l);
	return sli;
}

// slice a4(i,j,:,:)
template <class T, SLS_IF(is_scalar<T>())>
void slice34(Jcmat<T> &sli, const Cmat4d<T> &a4, Long_I i, Long_I j)
{
	Long N1 = a4.n1(), N2 = a4.n2(), N3 = a4.n3(), N4 = a4.n4();
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= N1 || j < 0 || j >= N2)
		SLS_ERR("out of bound!");
#endif
	Long N1N2 = N1 * N2;
	sli.set(a4.ptr() + i + N1*j, N3, N4, N1N2, N1N2*N3);
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat<T> slice34(const Cmat4d<T> &a4, Long_I i, Long_I j)
{
	Jcmat<T> sli;
	slice34(sli, a4, i, j);
	return sli;
}

// slice Jmat3d from Cmat3d
template <class T, SLS_IF(is_scalar<T>())>
void slice(Jcmat3d<T> &sli, const Cmat3d<T> &a3,
	Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k, Long_I N3)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i + N1 > a3.n1() || j < 0 ||
		j + N2 > a3.n2() || k < 0 || k + N2 > a3.n3())
		SLS_ERR("out of bound!");
#endif
	sli.set(&a3(i, j, k), N1, N2, N3, 1, a3.n1(), a3.n1()*a3.n2());
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat3d<T> slice(const Cmat3d<T> &a3,
	Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k, Long_I N3)
{
	Jcmat3d<T> sli;
	slice(sli, a3, i, N1, j, N2, k, N3);
	return sli;
}

// slice Jcmat4d from Cmat4d
template <class T, SLS_IF(is_scalar<T>())>
void slice(Jcmat4d<T> &sli, const Cmat4d<T> &a4,
	Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k, Long_I N3, Long_I l, Long_I N4)
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i + N1 > a4.n1() || j < 0 || j + N2 > a4.n2() ||
		k < 0 || k + N3 > a4.n3() || l < 0 || l + N4 > a4.n4())
		SLS_ERR("out of bound!");
#endif
	Long n1n2 = a4.n1() * a4.n2();
	sli.set(&a4(i, j, k, l), N1, N2, N3, N4, 1, a4.n1(), n1n2, n1n2*a4.n3());
}

template <class T, SLS_IF(is_scalar<T>())>
Jcmat4d<T> slice(const Cmat4d<T> &a4,
	Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k, Long_I N3, Long_I l, Long_I N4)
{
	Jcmat4d<T> sli;
	slice(sli, a4, i, N1, j, N2, k, N3, l, N4);
	return sli;
}

} // namespace slisc
