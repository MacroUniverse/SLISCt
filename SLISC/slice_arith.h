// slicing arithmetics
// slice() is for same dimensional container
// slice1() is for 1st dimension from heigher dimensional container
// slice12() is for first 2 dimensions from heigher dimensional container
// etc...

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
void slice_vec(Svector<T> &sli, Tdens &a, Long_I start, Long_I n)
{
#ifdef SLS_CHECKBOUNDS
    if (start < 0 || start + n > a.size())
        SLS_ERR("out of bound!");
#endif
    sli.set(a.ptr() + start, n);
}

template <class Tdens, class T = contain_type<Tdens>,
    SLS_IF(is_dense<Tdens>())>
Svector<T> slice_vec(Tdens &a, Long_I start, Long_I n)
{
    Svector<T> sli;
    slice_vec(sli, a, start, n);
    return sli;
}

template <class Tdens, class T = contain_type<Tdens>,
    SLS_IF(is_dense<Tdens>())>
inline void slice_vec(Svector_c<T> &sli, const Tdens &a, Long_I start, Long_I n)
{
#ifdef SLS_CHECKBOUNDS
    if (start < 0 || start + n > a.size())
        SLS_ERR("out of bound!");
#endif
    sli.set(a.ptr() + start, n);
}

template <class Tdens, class T = contain_type<Tdens>,
    SLS_IF(is_dense<Tdens>())>
inline Svector_c<T> slice_vec(const Tdens &a, Long_I start, Long_I n)
{
    Svector_c<T> sli;
    slice_vec(sli, a, start, n);
    return sli;
}

// slice a Dvector<> from any dense container using single index
template <class Tdens, class T = contain_type<Tdens>,
    SLS_IF(is_dense<Tdens>())>
void slice_vec(Dvector<T> &sli, const Tdens &a, Long_I start, Long_I n, Long_I step)
{
#ifdef SLS_CHECKBOUNDS
    if (start < 0 || start + n * step > a.size())
        SLS_ERR("out of bound!");
#endif
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
void slice2(Svector<T> &sli, Tmat &a, Long_I row)
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
Svector<T> slice2(Tmat &a, Long_I row)
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

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice2(Dvector<T> &sli, const Tmat &a, Long_I row, Long_I col, Long_I Ncol)
{
    Long Nr = a.n1();
#ifdef SLS_CHECK_BOUNDS
    if (row < 0 || row >= Nr)
        SLS_ERR("out of bound!");
    if (col + Ncol > a.n2())
        SLS_ERR("out of bound!");
#endif
    sli.set(a.ptr() + row + Nr*col, Ncol, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
Dvector<T> slice2(const Tmat &a, Long_I row, Long_I col, Long_I Ncol)
{
    Dvector<T> sli;
    slice2(sli, a, row, col, Ncol);
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

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline void slice2(Dvector<T> &sli, const Tmat &a, Long_I row,
        Long_I col, Long_I N2, Long_I k = 0)
{
    Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
    if (row < 0 || row >= Nr || col < 0 || col + N2 > Nc || N2 < 0)
        SLS_ERR("out of bound!");
#endif
    sli.set(&a(row, col, k), Nc, Nr);
}

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline Dvector<T> slice2(const Tmat &a, Long_I row,
        Long_I col, Long_I N2, Long_I k = 0)
{
    Dvector<T> sli;
    slice2(sli, a, row, col, N2, k);
    return sli;
}

// slice a col from a col-major matrix
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice1(Svector<T> &sli, Tmat &a, Long_I col)
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
Svector<T> slice1(Tmat &a, Long_I col)
{
    Svector<T> sli;
    slice1(sli, a, col);
    return sli;
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
void slice1(Svector_c<T> &sli, const Tmat &a, Long_I col)
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
Svector_c<T> slice1(const Tmat &a, Long_I col)
{
    Svector_c<T> sli;
    slice1(sli, a, col);
    return sli;
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
    inline void slice1(Svector<T> &sli, Tmat &a, Long_I row, Long_I N1, Long_I col)
{
    Long Nr = a.n1();
#ifdef SLS_CHECK_BOUNDS
    if (col < 0 || col >= a.n2())
        SLS_ERR("out of bound!");
    if (N1 > Nr)
        SLS_ERR("out of bound!");
#endif
    sli.set(a.ptr() + col * Nr + row, N1);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat<Tmat>() && is_cmajor<Tmat>())>
    inline Svector<T> slice1(Tmat &a, Long_I row, Long_I N1, Long_I col)
{
    Svector<T> sli;
    slice1(sli, a, row, N1, col);
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
void slice1(Svector<T> &sli, Tmat &a, Long_I col, Long_I k = 0)
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
inline Svector<T> slice1(Tmat &a, Long_I col, Long_I k)
{
    Svector<T> sli;
    slice1(sli, a, col, k);
    return sli;
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline void slice1(Svector_c<T> &sli, const Tmat &a, Long_I col, Long_I k = 0)
{
    Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_BOUNDS
    if (col < 0 || col >= Nc)
        SLS_ERR("out of bound!");
#endif
    sli.set(a.ptr() + Nr * Nc*k + Nr * col, Nr);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline Svector_c<T> slice1(const Tmat &a, Long_I col, Long_I k)
{
    Svector<T> sli;
    slice1(sli, a, col, k);
    return sli;
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline void slice1(Svector<T> &sli, Tmat &a,
    Long_I row, Long_I N1, Long_I col, Long_I k = 0)
{
#ifdef SLS_CHECK_BOUNDS
    if (N1 < 0 || row < 0 || row + N1 > a.n1())
        SLS_ERR("out of bound!");
    if (col < 0 || col >= a.n2())
        SLS_ERR("out of bound!");
    if (k < 0 || k >= a.n3())
        SLS_ERR("out of bound!");
#endif
    sli.set(&a(row, col, k), N1);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat3<Tmat>() && is_cmajor<Tmat>())>
inline Svector<T> slice1(Tmat &a,
    Long_I row, Long_I N1, Long_I col, Long_I k = 0)
{
    Svector<T> sli;
    slice1(sli, a, row, N1, col, k);
    return sli;
}

// slice a col from a col-major 4d array
template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat4<Tmat>() && is_cmajor<Tmat>())>
inline void slice1(Svector<T> &sli, Tmat &a, Long_I j, Long_I k = 0, Long_I l = 0)
{
    Long N1 = a.n1(), N2 = a.n2(), N3 = a.n3(), N4 = a.n4();
#ifdef SLS_CHECK_BOUNDS
    if (j < 0 || j >= N2 || k < 0 || k >= N3 || l < 0 || l >= N4)
        SLS_ERR("out of bound!");
#endif
    Long N1N2 = N1 * N2, N1N2N3 = N1N2 * N3;
    sli.set(a.ptr() + N1 * j + N1N2 * k + N1N2N3 * l, N1);
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense_mat4<Tmat>() && is_cmajor<Tmat>())>
inline Svector<T> slice1(Tmat &a, Long_I j, Long_I k = 0, Long_I l = 0)
{
    Svector<T> sli;
    slice1(sli, a, j, k, l);
    return sli;
}

// slice a Dmat from a dense mat
// only works for column major for now

template <class T, SLS_IF(is_scalar<T>())>
void slice(Scmat_c<T> &sli, const Cmat<T> &a, Long_I j, Long_I N2)
{
    sli.set(&a(0, j), a.n1(), N2);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat_c<T> slice(const Cmat<T> &a, Long_I j, Long_I N2)
{
    Scmat_c<T> sli;
    slice(sli, a, j, N2);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice(Scmat<T> &sli, Cmat<T> &a, Long_I j, Long_I N2)
{
    sli.set(&a(0, j), a.n1(), N2);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice(Cmat<T> &a, Long_I j, Long_I N2)
{
    Scmat<T> sli;
    slice(sli, a, j, N2);
    return sli;
}

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_cmajor<Tmat>() && is_dense_mat<Tmat>())>
void slice(Dcmat<T> &sli, const Tmat &a,
    Long_I i, Long_I Nr, Long_I j, Long_I Nc)
{
    sli.set(&a(i, j), Nr, Nc, a.n1());
}

template <class Tmat, class T = contain_type<Tmat>,
    SLS_IF(is_cmajor<Tmat>() && is_dense_mat<Tmat>())>
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
Dcmat<T> slice(const Dcmat<T> &a, Long_I i,
    Long_I Nr, Long_I j, Long_I Nc)
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
void slice12(Scmat_c<T> &sli, const Cmat3d<T> &a3, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
    if (k < 0 || k >= a3.n3())
        SLS_ERR("out of bound!");
#endif
    sli.set_size(a3.n1(), a3.n2());
    sli.set_ptr(a3.ptr() + sli.size() * k);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat_c<T> slice12(const Cmat3d<T> &a3, Long_I k)
{
    Scmat_c<T> sli;
    slice12(sli, a3, k);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat<T> &sli, Cmat3d<T> &a3, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
    if (k < 0 || k >= a3.n3())
        SLS_ERR("out of bound!");
#endif
    sli.set_size(a3.n1(), a3.n2());
    sli.set_ptr(a3.ptr() + sli.size() * k);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice12(Cmat3d<T> &a3, Long_I k)
{
    Scmat<T> sli;
    slice12(sli, a3, k);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat_c<T> &sli, const Cmat3d<T> &a3,
    Long_I j, Long_I N2, Long_I k)
{
    Long N1 = a3.n1();
    Long N1N2 = N1 * a3.n2();
    sli.set(a3.ptr() + N1 * j + N1N2 * k, N1, N2);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat_c<T> slice12(const Cmat3d<T> &a3,
    Long_I j, Long_I N2, Long_I k)
{
    Scmat_c<T> sli;
    slice12(sli, a3, j, N2, k);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat<T> &sli, Cmat3d<T> &a3, Long_I j, Long_I N2, Long_I k)
{
    Long N1 = a3.n1();
    Long N1N2 = N1 * a3.n2();
    sli.set(a3.ptr() + N1 * j + N1N2 * k, N1, N2);
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice12(Cmat3d<T> &a3,    Long_I j, Long_I N2, Long_I k)
{
    Scmat<T> sli;
    slice12(sli, a3, j, N2, k);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice12(Dcmat<T> &sli, const Cmat3d<T> &a3,
    Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k)
{
    sli.set(&a3(i, j, k), N1, N2, a3.n1());
}

template <class T, SLS_IF(is_scalar<T>())>
Dcmat<T> slice12(const Cmat3d<T> &a3,
    Long_I i, Long_I N1, Long_I j, Long_I N2, Long_I k)
{
    Dcmat<T> sli;
    slice12(sli, a3, i, N1, j, N2, k);
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
void slice12(Scmat_c<T> &a, const Cmat4d<T> &a4, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
    if (k < 0 || k >= a4.n3())
        SLS_ERR("out of bound!");
#endif
    a.set_size(a4.n1(), a4.n2());
    a.set_ptr(a4.ptr() + a.size() * (a4.n3()*l  + k));
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat_c<T> slice12(const Cmat4d<T> &a4, Long_I k, Long_I l)
{
    Scmat_c<T> sli;
    slice12(sli, a4, k, l);
    return sli;
}

template <class T, SLS_IF(is_scalar<T>())>
void slice12(Scmat<T> &a, Cmat4d<T> &a4, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
    if (k < 0 || k >= a4.n3())
        SLS_ERR("out of bound!");
#endif
    a.set_size(a4.n1(), a4.n2());
    a.set_ptr(a4.ptr() + a.size() * (a4.n3()*l + k));
}

template <class T, SLS_IF(is_scalar<T>())>
Scmat<T> slice12(Cmat4d<T> &a4, Long_I k, Long_I l)
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
        j + N2 > a3.n2() || k < 0 || k + N3 > a3.n3())
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
