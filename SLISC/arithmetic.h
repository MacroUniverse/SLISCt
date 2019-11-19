#pragma once
#include "global.h"
#include "meta.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "cmat3d.h"
#include "fixsize.h"
#include "ptr_arith.h"

namespace slisc {

// === get vec/mat properties ===

// check if vec/mat sizes are the same
template <class T1, class T2, SLS_IF(
    is_contain<T1>() && is_contain<T2>() &&
	ndims<T1>() == 1 && ndims<T2>() == 1)>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
    return v1.size() == v2.size();
}

template <class T1, class T2, SLS_IF(
    is_contain<T1>() && is_contain<T2>() &&
	ndims<T1>() == 2 && ndims<T2>() == 2)>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
    return v1.n1() == v2.n1() && v1.n2() == v2.n2();
}

template <class T1, class T2, SLS_IF(
    is_contain<T1>() && is_contain<T2>() &&
	ndims<T1>() == 3 && ndims<T2>() == 3)>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
        return v1.n1() == v2.n1() && v1.n2() == v2.n2()
            && v1.n3() == v2.n3();
}

template <class T1, class T2, SLS_IF(
    is_contain<T1>() && is_contain<T2>() &&
	ndims<T1>() != ndims<T2>())>
Bool shape_cmp(const T1 &v1, const T2 &v2)
{
        return false;
}

// operator== for slisc containers
template <class T1, class T2, SLS_IF(
    is_dense<T1>() && is_dense<T2>() &&
    is_same_major<T1, T2>())>
Bool operator==(const T1 &v1, const T2 &v2)
{
    return shape_cmp(v1, v2) &&
        equals_to_vv(v1.ptr(), v2.ptr(), v2.size());
}

// for Mcoo and Mcooh matrices
template <class T1, class T2, SLS_IF(
    (is_MatCoo<T1>() || is_MatCooH<T1>()) &&
     is_same_contain<T1, T2>())>
Bool operator==(const T1 &v1, const T2 &v2)
{
    Long Nnz = v1.nnz();
    if (!shape_cmp(v1, v2) || Nnz != v2.nnz())
        return false;
    T1 u1(v1.n1(), v1.n2(), Nnz);
    u1 = v1;
    T2 u2(v2.n1(), v2.n2(), Nnz);
    u2 = v2;
    u1.sort_r(); u2.sort_r();
    for (Long i = 0; i < Nnz; ++i) {
        if (u1[i] != u2[i] || u1.row(i) != u2.row(i) ||
            u1.col(i) != u2.col(i))
            return false;
    }
    return true;
}

template <class T1, class T2, SLS_IF(
    is_Dvector<T1>() && is_Dvector<T2>())>
Bool operator==(const T1 &v1, const T2 &v2)
{
    if (!shape_cmp(v1, v2))
        return false;
    for (Long i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i])
            return false;
    }
    return true;
}

// for all other (2d) matrices not covered by the above definition
// might not be efficient
template <class T1, class T2, SLS_IF(
    ndims<T1>() == 2 && ndims<T2>() == 2 &&
    ! (is_dense<T1>() && is_dense<T2>() &&
        is_same_major<T1, T2>()) &&
    ! ((is_MatCoo<T1>() || is_MatCooH<T1>()) &&
        is_same_contain<T1, T2>()))>
Bool operator==(const T1 &v1, const T2 &v2)
{
    if (!shape_cmp(v1, v2))
        return false;
    for (Long i = 0; i < v1.n1(); ++i)
        for (Long j = 0; j < v1.n2(); ++j)
            if (v1(i, j) != v2(i, j))
                return false;
    return true;
}

// all other 3d matrices
template <class T1, class T2, SLS_IF(
    ndims<T1>() == 3 && ndims<T2>() == 3 &&
    !(is_dense<T1>() && is_dense<T2>() &&
    is_same_major<T1, T2>()))>
inline Bool operator==(const T1 &v1, const T2 &v2)
{
    if (!shape_cmp(v1, v2))
        return false;
    for (Long k = 0; k < v1.n3(); ++k)
        for (Long j = 0; j < v1.n2(); ++j)
            for (Long i = 0; i < v1.n1(); ++i)
                if (v1(i, j, k) != v2(i, j, k))
                    return false;
    return true;
}

template <class Tv, class Ts, SLS_IF(is_dense<Tv>() && is_scalar<Ts>())>
inline Bool operator==(const Tv &v, const Ts &s)
{
    return equals_to_vs(v.ptr(), s, v.size());
}

template <class Tv, class Ts, SLS_IF(is_dense<Tv>() && !is_contain<Ts>())>
inline Bool operator==(const Ts &s, const Tv &v)
{ return v == s; }

template <class Tv, class Ts, SLS_IF(
    ndims<Tv>() == 1 && !is_dense<Tv>() && is_scalar<Ts>())>
inline Bool operator==(const Tv &v, const Ts &s)
{
    for (Long i = 0; i < v.size(); ++i) {
        if (v[i] != s)
            return false;
    }
    return true;
}

// operator!= for slisc containers
template <class T1, class T2, SLS_IF(
    is_contain<T1>() && is_scalar<T2>() ||
    is_contain<T2>() && is_scalar<T1>() ||
    ndims<T1>() > 0 && ndims<T1>() == ndims<T2>())>
Bool operator!=(const T1 &v1, const T2 &v2)
{
    return !(v1 == v2);
}

// copy one row of dense matrix to a dense vector
template <class Tvec, class Tmat,
    SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_row(Tvec &v, const Tmat &a, Long_I row)
{
    Long Nc = a.n2();
#ifdef SLS_CHECK_SHAPE
    if (v.size() != Nc) {
        SLS_ERR("wrong shape!");
    }
#endif
    if (is_rmajor<Tmat>()) { // row major
        veccpy(v.ptr(), a.ptr() + Nc*row, Nc);
    }
    else if (is_cmajor<Tmat>()) { // column major
        Long Nr = a.n1();
        auto p = a.ptr() + row;
        for (Long i = 0; i < Nc; ++i) {
            v[i] = *p;
            p += Nr;
        }
    }
    else
        SLS_ERR("unknown!");
}

// copy a dense vector to one row of dense matrix
template <class Tvec, class Tmat,
    SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_row(Tmat &a, const Tvec &v, Long_I row)
{
    Long Nc = a.n2();
#ifdef SLS_CHECK_SHAPE
    if (v.size() != Nc) {
        SLS_ERR("wrong shape!");
    }
#endif
    if (is_rmajor<Tmat>()) { // row major
        veccpy(a.ptr() + Nc * row, v.ptr(), Nc);
    }
    else if (is_cmajor<Tmat>()) { // column major
        Long Nr = a.n1();
        auto p = a.ptr() + row;
        for (Long i = 0; i < Nc; ++i) {
            *p = v[i];
            p += Nr;
        }
    }
    else
        SLS_ERR("unknown!");
}

// copy one column of dense matrix to a dense vector
template <class Tvec, class Tmat,
    SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_col(Tvec &v, const Tmat &a, Long_I col)
{
#ifdef SLS_CHECK_SHAPE
    if (v.size() != a.n1()) {
        SLS_ERR("wrong shape!");
    }
#endif
    Long Nr = a.n1();
    if (is_cmajor<Tmat>()) { // column major
        veccpy(v.ptr(), a.ptr() + Nr*col, Nr);
    }
    else if (is_rmajor<Tmat>()) { // row major
        Long Nc = a.n2();
        auto p = a.ptr() + col;
        for (Long i = 0; i < Nr; ++i) {
            v[i] = *p;
            p += Nc;
        }
    }
    else
        SLS_ERR("unknown!");
}

// copy a dense vector to one column of dense matrix
template <class Tvec, class Tmat,
    SLS_IF(is_dense_vec<Tvec>() && is_dense_mat<Tmat>())>
inline void copy_col(Tmat &a, const Tvec &v, Long_I col)
{
#ifdef SLS_CHECK_SHAPE
    if (v.size() != a.n1()) {
        SLS_ERR("wrong shape!");
    }
#endif
    Long Nr = a.n1();
    if (is_cmajor<Tmat>()) { // column major
        veccpy(a.ptr() + Nr * col, v.ptr(), Nr);
    }
    else if (is_rmajor<Tmat>()) { // row major
        Long Nc = a.n2();
        auto p = a.ptr() + col;
        for (Long i = 0; i < Nr; ++i) {
            *p = v[i];
            p += Nc;
        }
    }
    else
        SLS_ERR("unknown!");
}

// vector/matrix stride (step2 is leading demension for matrix)
template <class T, SLS_IF(is_dense_vec<T>())>
Long step1(const T &v) {
	return 1;
}

template <class T, SLS_IF(is_Dvector<T>())>
Long step1(const T &v) {
	return v.step();
}

template <class T, SLS_IF(is_dense<T>() && is_cmajor<T>())>
Long step1(const T &v) {
	return 1;
}

template <class T, SLS_IF(is_dense<T>() && is_cmajor<T>())>
Long step2(const T &v) {
	return v.n1();
}

template <class T, SLS_IF(is_Dcmat<T>())>
Long step2(const T &v) {
	return v.lda();
}

// sum of container elements
template <class T, SLS_IF(is_dense<T>())>
inline const auto sum(const T &v)
{
    return sum_v(v.ptr(), v.size());
}

// sum of abs of elements
template <class T, SLS_IF(is_dense<T>())>
inline const auto sum_abs(const T &v)
{
    return sum_abs_v(v.ptr(), v.size());
}

template <class T, SLS_IF(is_dense<T>())>
inline const auto max(const T &v)
{
    return max_v(v.ptr(), v.size());
}

template <class T, SLS_IF(is_Dvector<T>())>
inline const auto max(const T &v)
{
    return max_v(v.ptr(), v.size(), v.step());
}

// return max(abs(a(:))
template <class T, SLS_IF(is_dense<T>())>
inline const auto max_abs(const T &v)
{
    return max_abs_v(v.ptr(), v.size());
}

template <class T, SLS_IF(is_Dvector<T>())>
inline const auto max_abs(const T &v)
{
    return max_abs_v(v.ptr(), v.size(), v.step());
}

// return min(abs(a(:))
template <class T, SLS_IF(is_dense<T>())>
inline const auto min_abs(const T &v)
{
    return min_abs_v(v.ptr(), v.size());
}

template <class T, SLS_IF(is_Dvector<T>())>
inline const auto min_abs(const T &v)
{
    return min_abs_v(v.ptr(), v.size(), v.step());
}

template <class T, SLS_IF(is_dense<T>())>
inline const auto max(Long_O ind, const T &v)
{
    Long i, N{ v.size() };
    auto val = v[0];
    for (i = 1; i < N; ++i) {
        if (val < v[i]) {
            val = v[i]; ind = i;
        }
    }
    return val;
}

// s = norm2(v)  (|v1|^2 + |v2|^2 + ...)
template <class T, SLS_IF(is_dense<T>() || is_Dvector<T>())>
const auto norm2(const T &v)
{
    Long N = v.size();
    auto s2 = ABS2(v[0]);
    for (Long i = 1; i < N; ++i)
        s2 += ABS2(v[i]);
    return s2;
}

template <class T, SLS_IF(is_Dcmat<T>())>
const auto norm2(const T &a)
{
    auto p = a.ptr();
    Long Nr = a.n1(), lda = a.lda();
    rm_comp<contain_type<T>> s2 = 0;
    for (Long j = 0; j < a.n2(); ++j) {
        for (Long i = 0; i < Nr; ++i)
            s2 += ABS2(p[i]);
        p += lda;
    }
    return s2;
}

template <class T>
const auto norm(const T &v, SLS_IF(is_dense<T>()))
{
    return sqrt(norm2(v));
}

// === matrix manipulation ===

// does not work for integers

template <class Tv, class T1, class T2, SLS_IF(is_dense<Tv>())>
inline void linspace(Tv &v, const T1 &first, const T2 &last)
{
    linspace_vss(v.ptr(), (contain_type<Tv>)first, (contain_type<Tv>)last, v.size());
}

template <class Tv, SLS_IF(ndims<Tv>() == 1)>
inline void reorder(Tv &v, VecLong_I order)
{
    Long N = v.size();
#ifdef SLS_CHECK_SHAPE
    if (order.size() != N)
        SLS_ERR("wrong shape!");
#endif
    Vector<contain_type<Tv>> u(N);
    for (Long i = 0; i < N; ++i) {
        u[i] = v[order[i]];
    }
    for (Long i = 0; i < N; ++i) {
        v[i] = u[i];
    }
}

// === vectorized math functions ===

template <class T, SLS_IF(is_dense<T>())>
inline void sqrt(T &v)
{
    sqrt_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_dense<T1>())>
inline void sqrt(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    sqrt_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class Ts, SLS_IF(
    is_dense<T>() && is_scalar<Ts>())>
inline void pow(T &v, const Ts &s)
{
    pow_vs(v.ptr(), s, v.size());
}

template <class T, class T1, class T2, SLS_IF(
    is_dense<T>() && is_same_contain<T, T1>() && is_scalar<T2>())>
inline void pow(T &v, const T1 &v1, const T2 &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    pow_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void sin(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    sin_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void cos(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    cos_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, SLS_IF(is_fpt_dense<T>())>
inline void exp(T &v)
{
    exp_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void exp(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    exp_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1, SLS_IF(is_dense<T>() && is_same_contain<T, T1>())>
inline void tan(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    tan_vv(v.ptr(), v1.ptr(), v1.size());
}

// === vector/matrix arithmetics ===

// flip(v)
template <class T, SLS_IF(is_dense_vec<T>())>
inline void flip(T &v)
{
    flip_v(v.ptr(), v.size());
}

// v = flip(v)
template <class T, class T1, SLS_IF(
    is_dense_vec<T>() && is_dense_vec<T1>()
)>
inline void flip(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    flip(v.ptr(), v1.ptr(), v1.size());
}

// matrix transpose

// trans(v)
template <class T, SLS_IF(is_dense_mat<T>())>
inline void trans(T &v)
{
#ifdef SLS_CHECK_SHAPE
    if (v.n1() != v.n2())
        SLS_ERR("illegal shape!");
#endif
    for (Long i = 0; i < v.n1(); ++i)
        for (Long j = 0; j < i; ++j)
            swap(v(i, j), v(j, i));
}

// v = trans(v)
// TODO: this is inefficient
template <class T, class T1, SLS_IF(
    is_dense_mat<T>() && is_dense_mat<T1>()
)>
inline void trans(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (v.n1() != v1.n2() || v.n2() != v1.n1())
        SLS_ERR("wrong size!");
#endif
    for (Long i = 0; i < v.n1(); ++i)
        for (Long j = 0; j < v.n2(); ++j)
            v(i, j) = v1(j, i);
}

// hermitian conjugate

// her(v)
template <class T, SLS_IF(is_comp_dense<T>() && is_dense_mat<T>())>
inline void her(T &v)
{
#ifdef SLS_CHECK_SHAPE
    if (v.n1() != v.n2()) SLS_ERR("illegal shape!");
#endif
    // diagonal
    for (Long i = 0; i < v.n1(); ++i)
        v(i, i) = conj(v(i, i));
    // off-diagonal
    for (Long i = 0; i < v.n1(); ++i) {
        for (Long j = 0; j < i; ++j) {
            contain_type<T> temp = v(i, j);
            v(i, j) = conj(v(j, i));
            v(j, i) = conj(temp);
        }
    }
}

// v = her(v)
template <class T, class T1, SLS_IF(
    is_dense_mat<T>() && is_dense_mat<T1>() &&
    is_comp_dense<T>() && is_comp_dense<T1>() &&
    type_num<contain_type<T>>() >= type_num<contain_type<T1>>()
)>
inline void her(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (v.n1() != v1.n2() || v.n2() != v1.n1())
        SLS_ERR("wrong shape!");
#endif
    for (Long i = 0; i < v.n1(); ++i) {
        for (Long j = 0; j < v.n2(); ++j)
            v(i, j) = conj(v1(j, i));
    }
}

// v += v

template <class T, class T1, SLS_IF(is_dense<T>() && is_dense<T1>())>
inline void operator+=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    plus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v -= v

template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void operator-=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    minus_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

template <class T, class T1, SLS_IF(
is_contain<T>() && is_contain<T1>() &&
(!is_dense<T>() || !is_dense<T1>()) &&
ndims<T>() == 1 && ndims<T1>() == 1 &&
is_promo<contain_type<T>, contain_type<T1>>())>
inline void operator-=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    for (Long i = 0; i < v.size(); ++i)
        v[i] -= v1[i];
}

// v *= v

template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void operator*=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    times_equals_vv(v.ptr(), v1.ptr(), v1.size());
}



// v /= v

template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void operator/=(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    divide_equals_vv(v.ptr(), v1.ptr(), v1.size());
}

// v += s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator+=(T &v, const Ts &s)
{
    plus_equals_vs(v.ptr(), s, v.size());
}

// v -= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator-=(T &v, const Ts &s)
{
    minus_equals_vs(v.ptr(), s, v.size());
}

// v *= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator*=(T &v, const Ts &s)
{
    times_equals_vs(v.ptr(), s, v.size());
}

template <class T, class Ts, SLS_IF(is_Dvector<T>() && is_scalar<Ts>())>
inline void operator*=(T &v, const Ts &s)
{
    times_equals_vs(v.ptr(), s, v.size(), v.step());
}

template <class T, class T1, SLS_IF(
    is_scalar<T>() && is_promo<T, T1>())>
inline void operator*=(Dcmat<T> &v, Dcmat<T1> &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    for (Long j = 0; j < v.n2(); ++j)
        times_equals_vv(&v(0, j), &v1(0, j), v.n1());
}

// v /= s

template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void operator/=(T &v, const Ts &s)
{
    divide_equals_vs(v.ptr(), s, v.size());
}

template <class T, class Ts, SLS_IF(is_Dvector<T>() && is_scalar<Ts>())>
inline void operator/=(T &v, const Ts &s)
{
    divide_equals_vs(v.ptr(), s, v.size(), v.step());
}

// v %= s
template <class T, class Ts, SLS_IF(is_dense<T>() && is_scalar<Ts>())>
inline void rem(T &v, const Ts &s)
{
    rem_vs(v.ptr(), s, v.size());
}

// v = v % s

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void rem(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    rem_vvs(v.ptr(), v1.ptr(), s, v.size());
}

// mod(v, s)
template <class T, class T1, SLS_IF(
    is_dense<T>() && is_scalar<T1>())>
inline void mod(T &v, const T1 &s)
{
    mod_vs(v.ptr(), s, v.size());
}

// v = mod(v, s)
template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void mod(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    mod_vvs(v.ptr(), v1.ptr(), s, v.size());
}

// TODO : mod(v, s, v1)
// TODO : mod(v, v1, v2)

// v = v + s

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>()
)>
inline void Plus(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Plus(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    plus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v + v

template <class T, class T1, class T2, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_dense<T2>()
)>
inline void Plus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v1, v2) || !shape_cmp(v, v1)) SLS_ERR("wrong shape!");
#endif
    plus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// -v inplace
template <class T, SLS_IF(is_dense<T>())>
inline void Minus(T &v)
{
    minus_v(v.ptr(), v.size());
}

// v = -v

template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void Minus(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    minus_vv(v.ptr(), v1.ptr(), v1.size());
}

// v = s - v

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Minus(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    minus_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v - s

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Minus(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    minus_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v - v

template <class T, class T1, class T2, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_dense<T2>())>
inline void Minus(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v1, v2) || !shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    minus_vvv(v.ptr(), v1.ptr(), v2.ptr(), v1.size());
}

// v = v * s

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Times(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s * v

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Times(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    times_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = v * v

template <class T, class T1, class T2, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_dense<T2>() &&
    is_same_major<T, T1>() && is_same_major<T, T2>())>
inline void Times(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v1, v2) || !shape_cmp(v, v1)) SLS_ERR("wrong shape!");
#endif
    times_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// v = v / s

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Divide(T &v, const T1 &v1, const Ts &s)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    divide_vvs(v.ptr(), v1.ptr(), s, v1.size());
}

// v = s / v

template <class T, class T1, class Ts, SLS_IF(
    is_dense<T>() && is_dense<T1>() &&
    is_same_major<T, T1>() && is_scalar<Ts>())>
inline void Divide(T &v, const Ts &s, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    divide_vsv(v.ptr(), s, v1.ptr(), v1.size());
}

// v = v / v

template <class T, class T1, class T2, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_dense<T2>() &&
    is_same_major<T, T1>() && is_same_major<T, T2>())>
inline void Divide(T &v, const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v1, v2) || !shape_cmp(v, v1))
        SLS_ERR("wrong shape!");
#endif
    divide_vvv(v.ptr(), v1.ptr(), v2.ptr(), v2.size());
}

// real(v)

template <class T, SLS_IF(is_comp_dense<T>())>
inline void real(T &v)
{
    real_v(v.ptr(), v.size());
}

// v = real(v)

template <class T, class T1,
    SLS_IF(is_dense<T>() && is_comp_dense<T1>() &&
        is_same_major<T, T1>())>
inline void real(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    real_vv(v.ptr(), v1.ptr(), v1.size());
}

// imag(v)
template <class T, SLS_IF(is_comp_dense<T>())>
inline void imag(T &v)
{
    imag_v(v.ptr(), v.size());
}

// v = imag(v)
template <class T, class T1,
    SLS_IF(is_dense<T>() && is_comp_dense<T1>() &&
        is_same_major<T, T1>())>
inline void imag(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    imag_vv(v.ptr(), v1.ptr(), v1.size());
}

// abs(v)
template <class T, SLS_IF(is_dense<T>())>
inline void abs(T &v)
{
    abs_v(v.ptr(), v.size());
}

// v = abs(v)
template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void abs(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    abs_vv(v.ptr(), v1.ptr(), v1.size());
}

// v = abs(v)^2
template <class T, class T1, SLS_IF(
    is_dense<T>() && is_dense<T1>() && is_same_major<T, T1>())>
inline void abs2(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    abs2_vv(v.ptr(), v1.ptr(), v1.size());
}

// to_comp() deleted, use operator= instead

// conj(v)
template <class T, SLS_IF(is_comp_dense<T>())>
inline void conj(T &v)
{
    conj_v(v.ptr(), v.size());
}

template <class T, class T1, SLS_IF(
    is_comp_dense<T>() && is_comp_dense<T1>() &&
    is_same_major<T, T1>())>
inline void conj(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("wrong size!");
#endif
    conj_vv(v.ptr(), v1.ptr(), v1.size());
}

// dot products ( sum conj(v1[i])*v2[i] )

// s = dot(v, v)
template <class T1, class T2, SLS_IF(is_dense<T1>() && is_dense<T2>())>
inline auto dot(const T1 &v1, const T2 &v2)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v1, v2))
        SLS_ERR("wrong shape!");
#endif
    return dot_vv(v1.ptr(), v2.ptr(), v2.size());
}

// outer product ( v(i,j) = conj(v1[i])*v2[j] )
// outprod(v, v, v)
template <class T, class T1, class T2, SLS_IF(
    is_dense_mat<T>() && is_dense_vec<T1>() && is_dense_vec<T2>())>
inline void outprod(T &v, const T1 &v1, const T2 &v2)
{
    SLS_ERR("TODO");
    /*Long i, j, N1{ v1.size() }, N2{ v2.size() };
    Comp *pc, v1_i;
    v.resize(N1, N2);
    for (i = 0; i < N1; ++i) {
        pc = v[i];
        v1_i = v1[i];
        for (j = 0; j < N2; ++j)
            pc[j] = v1_i*v2[j];
    }*/
}

template <class T, class T1, class T2, SLS_IF(
    is_dense_mat<T>() && is_dense_vec<T1>() && is_dense_vec<T2>()
)>
inline void outprod_par(T &v, const T1 &v1, const T2 &v2)
{
    SLS_ERR("TODO");
    /*Long N1{ v1.size() }, N2{ v2.size() };
    v.resize(N1, N2);
#pragma omp parallel for
    for (Long i = 0; i < N1; ++i) {
        Comp *pc, v1_i;
        pc = v.ptr(i);
        v1_i = v1[i];
        for (Long j = 0; j < N2; ++j)
            pc[j] = v1_i*v2[j];
    }*/
}

// matrix-vector multiplication
template <class T, class T1, class T2, SLS_IF(
    ndims<T>() == 1 && ndims<T2>() == 1 &&
    (is_dense_mat<T1>() || is_Dcmat<T1>()))>
inline void mul(T &y, const T1 &a, const T2 &x)
{
    Long Nr = a.n1(), Nc = a.n2();
#ifdef SLS_CHECK_SHAPE
    if (Nc != x.size() || y.size() != Nr)
        SLS_ERR("illegal shape!");
#endif
    for (Long i = 0; i < Nr; ++i) {
        y[i] = 0;
        for (Long j = 0; j < Nc; ++j)
            y[i] += a(i, j) * x[j];
    }
}

// matrix-vector multiplication with symmetric Doub matrix and Comp vectors (use MKL)
template <class T, class T1, class T2, SLS_IF(
    is_dense_vec<T>() && is_dense_mat<T1>() && is_dense_vec<T2>() &&
    is_Comp<contain_type<T>>() && is_Doub<contain_type<T1>>() &&
    is_Comp<contain_type<T2>>()
)>
inline void mul_sym(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
    if (a.n1() != a.n2() || x.size() != y.size() || x.size() != a.n1())
        SLS_ERR("wrong shape!");
#endif
#ifdef SLS_USE_CBLAS
    // do real part
    Long N = x.size();
    cblas_dsymv(CblasColMajor, CblasUpper, N, 1, a.ptr(),
        N, (Doub*)x.ptr(), 2, 0, (Doub*)y.ptr(), 2);
    // do imag part
    cblas_dsymv(CblasColMajor, CblasUpper, N, 1, a.ptr(),
        N, (Doub*)x.ptr() + 1, 2, 0, (Doub*)y.ptr() + 1, 2);
#else
    mul(y, a, x);
#endif
}

// matrix-vector multiplication with symmetric Dcmat and vectors of Doub and Comp (use MKL)
template <class T, class T1, class T2, class Ts = contain_type<T>,
    class Ts1 = contain_type<T1>, class Ts2 = contain_type<T2>, SLS_IF(
    is_dense_vec<T>() && is_Dcmat<T1>() && is_dense_vec<T2>() &&
    is_Doub<Ts>() && is_Doub<Ts1> && is_Doub<Ts1>
)>
inline void mul_sym(T &y, const T1 &a, const T2 &x)
{
    Long N = x.size();
#ifdef SLS_CHECK_SHAPE
    if (a.n1() != a.n2() || N != y.size() || N != a.n1())
        SLS_ERR("wrong shape!");
#endif
#ifdef SLS_USE_CBLAS
    if (is_Doub<Ts>()) {
        cblas_dsymv(CblasColMajor, CblasUpper, N, 1, a.ptr(),
            a.lda(), (Doub*)x.ptr(), 2, 0, (Doub*)y.ptr(), 2);
    }
	else
		SLS_ERR("not implemented");
#else
    mul(y, a, x);
#endif
}

// matrix-vector multiplication using cBLAS
template <class T, class T1, class T2,
    class Ts = contain_type<T>, class Ts1 = contain_type<T1>,
    class Ts2 = contain_type<T2>, SLS_IF(
		is_cmajor<T1>() &&
        (is_dense_vec<T>() || is_Dvector<T>()) &&
        (is_dense_mat<T1>() || is_Dcmat<T1>()) &&
        (is_dense_vec<T2>() || is_Dvector<T2>()) &&
        ((is_Doub<Ts>() && is_Doub<Ts1>() && is_Doub<Ts2>()) ||
         (is_Comp<Ts>() && is_Comp<Ts1>() && is_Comp<Ts2>()) ||
         (is_Comp<Ts>() && is_Doub<Ts1>() && is_Comp<Ts>())))>
inline void mul_gen(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
    if (x.size() != a.n2() || y.size() != a.n1())
        SLS_ERR("wrong shape!");
#endif
#ifdef SLS_USE_CBLAS
	Long N1 = a.n1(), N2 = a.n2(), lda, incx, incy;
    incy =  step1(y);
	lda = step2(a);
	incx = step1(x);

    if (is_Doub<Ts>() && is_Doub<Ts1>())
        cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N2, 1, (Doub *)a.ptr(),
            lda, (Doub *)x.ptr(), incx, 0, (Doub *)y.ptr(), incy);
    else if (is_Comp<Ts>() && is_Comp<Ts1>()) {
        Comp alpha(1), beta(0);
        cblas_zgemv(CblasColMajor, CblasNoTrans, N1, N2, &alpha, (Comp *)a.ptr(),
            lda, (Comp *)x.ptr(), incx, &beta, (Comp *)y.ptr(), incy);
    }
    else if (is_Comp<Ts>() && is_Doub<Ts1>()) {
        // do real part
        cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N2, 1, (Doub*)a.ptr(),
            lda, (Doub*)x.ptr(), 2*incx, 0, (Doub*)y.ptr(), 2*incy);
        // do imag part
        cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N2, 1, (Doub*)a.ptr(),
            lda, (Doub*)x.ptr() + 1, 2*incx, 0, (Doub*)y.ptr() + 1, 2*incy);
    }
#else
    mul(y, a, x);
#endif
}

// parallel version
template <class T, class T1, class T2, SLS_IF(
    is_dense_vec<T>() && is_dense_mat<T1>() && is_dense_vec<T2>()
)>
inline void mul_par(T &y, const T1 &a, const T2 &x)
{
#ifdef SLS_CHECK_SHAPE
    if (a.n2() != x.size() || y.size() != a.n1())
        SLS_ERR("illegal shape!");
#endif
    Long Nr_a = a.n1(), Nc_a = a.n2();
    vecset(y.ptr(), contain_type<T>(), Nr_a);
#pragma omp parallel for
    for (Long i = 0; i < Nr_a; ++i) {
        for (Long j = 0; j < Nc_a; ++j)
            y[i] += a(i, j) * x[j];
    }
}

// vector-matrix multiplication (row vector assumed)
template <class T, class T1, class T2, SLS_IF(
    is_dense_vec<T>() && is_dense_vec<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &x, const T2 &a)
{
    Long Nr_a = a.n1(), Nc_a = a.n2();
#ifdef SLS_CHECK_SHAPE
    if (x.size() != a.n1() || y.size() != Nc_a)
        SLS_ERR("illegal shape!");
#endif
    vecset(y.ptr(), contain_type<T>(), Nc_a);
    for (Long j = 0; j < Nc_a; ++j) {
        for (Long i = 0; i < Nr_a; ++i)
            y[j] += x[i] * a(i, j);
    }
}

// parallel version
template <class T, class T1, class T2, SLS_IF(
    is_dense_vec<T>() && is_dense_vec<T1>() && is_dense_mat<T2>()
)>
inline void mul_par(T &y, const T1 &x, const T2 &a)
{
    Long Nr_a = a.n1(), Nc_a = a.n2();
#ifdef SLS_CHECK_SHAPE
    if (x.size() != a.n1() || y.size() != Nc_a)
        SLS_ERR("illegal shape!");
#endif
    vecset(y.ptr(), contain_type<T>(), Nc_a);
#pragma omp parallel for
    for (Long j = 0; j < Nc_a; ++j) {
        for (Long i = 0; i < Nr_a; ++i)
            y[j] += x[i] * a(i, j);
    }
}

// matrix-matrix multiplication

template <class T, class T1, class T2, SLS_IF(
    is_dense_mat<T>() && is_dense_mat<T1>() && is_dense_mat<T2>()
)>
inline void mul(T &y, const T1 &a, const T2 &x)
{
    Long Nr_a = a.n1(), Nc_a = a.n2(), Nc_x = x.n2();
#ifdef SLS_CHECK_SHAPE
    if (a.n2() != x.n1() || y.n1() != Nr_a || y.n2() != Nc_x)
        SLS_ERR("illegal shape!");
#endif
    vecset(y.ptr(), contain_type<T>(), Nr_a*Nc_x);
    for (Long i = 0; i < Nr_a; ++i) {
        for (Long j = 0; j < Nc_x; ++j) {
            for (Long k = 0; k < Nc_a; ++k)
                y(i, j) += a(i, k) * x(k, j);
        }
    }
}

// using mkl
template <class T, class T1, class T2, 
	class Ts = contain_type<T>, class Ts1 = contain_type<T1>,
	class Ts2 = contain_type<T2>, SLS_IF(
    is_dense_mat<T>() && is_dense_mat<T1>() && is_dense_mat<T2>())>
inline void mul_gen(T &y, const T1 &a, const T2 &x)
{
    Long Nr_a = a.n1(), Nc_a = a.n2(), Nc_x = x.n2();
#ifdef SLS_CHECK_SHAPE
    if (a.n2() != x.n1() || y.n1() != Nr_a || y.n2() != Nc_x)
        SLS_ERR("illegal shape!");
#endif
#ifdef SLS_USE_CBLAS
    if (is_Doub<Ts>() && is_Doub<Ts1>() && is_Doub<Ts2>())
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nr_a, Nc_x, Nc_a, 1, (Doub*)a.ptr(), Nr_a, (Doub*)x.ptr(), Nc_a, 0, (Doub*)y.ptr(), Nr_a);
    else if (is_Comp<Ts>() && is_Comp<Ts1>() && is_Comp<Ts2>()) {
        Comp alpha(1,0), beta(0,0);
        // this might cause memory read violation, don't know why!
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nr_a, Nc_x, Nc_a, &alpha, a.ptr(), Nr_a, x.ptr(), Nc_a, &beta, y.ptr(), Nr_a);
    }
    else
        SLS_WARN("not implemented with cBLAS, using slow version");
		mul(y, a, x);
#else
    mul(y, a, x);
#endif
}

// === numerical integration ===

// indefinite integral;
// use cumsum(y)*dx instead
template <class T, class T1, SLS_IF(is_dense_vec<T>() && is_dense_vec<T1>())>
inline void cumsum(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(v, v1))
        SLS_ERR("illegal shape!");
#endif
        cumsum_vv(v.ptr(), v1.ptr(), v1.size());
}

// === others ===

// get all unique rows from a matrix
template <class Tmat, class Tmat1, class T = contain_type<Tmat>,
    class T1 = contain_type<Tmat1>, SLS_IF(
    ndims<Tmat>() == 2 && ndims<Tmat1>() ==2 &&
    is_cmajor<Tmat>() && is_cmajor<Tmat1>() &&
    is_promo<T, T1>())>
void uniq_rows(Tmat &a, const Tmat1 &a1)
{
    Long k = 0;
    Dvector<T1> sli1;
    a.resize(a1.n1(), a1.n2());
    for (Long i = 0; i < a1.n1(); ++i) {
        // check repeat
        Bool repeat = false;
        slice2(sli1, a1, i);
        for (Long j = 0; j < k; ++j) {
            if (slice2(a, j) == sli1) {
                repeat = true; break;
            }
        }
        if (repeat)
            continue;
        slice2(a, k) = sli1;
        ++k;
    }
    a.resize_cpy(k, a1.n2());
}

} // namespace slisc
