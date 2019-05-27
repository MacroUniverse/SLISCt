// copy data from one container to another
// includes container shape checking
// all `container = container` should be implemented using copy

#include "meta.h"

namespace slisc {
// memory set and copy
template<class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void vecset(T *dest, const T1 &val, Long_I n)
{
	T val0 = (T)val;
	for (T *p = dest; p < dest + n; ++p)
		*p = val0;
}

template<class T, class T1, SLS_IF(is_promo<T, T1>())>
inline void vecset(T *dest, const T1 &val, Long_I n, Long_I step)
{
	T val0 = (T)val;
	for (T *p = dest; p < dest + n*step; p += step)
		*p = val0;
}

template<class T>
inline void veccpy(T *dest, const T *src, Long_I n)
{
	memcpy(dest, src, n * sizeof(T));
}

template<class T, class T1, SLS_IF(is_promo<T,T1>())>
inline void veccpy(T *dest, const T1 *src, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = src[i];
}

template <class T, class T1, SLS_IF(
	is_dense_vec<T>() && is_dense_vec<T1>() ||
	is_dense<T>() && is_dense<T1>() &&
	is_same_major<T, T1>())>
inline void copy(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	veccpy(v.ptr(), v1.ptr(), v.size());
}

// copy from row-major to col-major
template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>() &&
	is_cmajor<T>() && is_rmajor<T1>())>
inline void copy(T &a, const T1 &a1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(a, a1))
		SLS_ERR("wrong shape!");
#endif
	Long Nr = a.nrows(), Nc = a.ncols(), k1 = 0;
	for (Long i = 0; i < Nr; ++i) {
		Long k = i;
		for (Long j = 0; j < Nc; ++j) {
			a[k] = a1[k1];
			k += Nc; ++k1;
		}
	}
}

// copy from col-major to row-major
template <class T, class T1, SLS_IF(
	is_dense_mat<T>() && is_dense_mat<T1>() &&
	is_rmajor<T>() && is_cmajor<T1>())>
inline void copy(T &a, const T1 &a1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(a, a1))
		SLS_ERR("wrong shape!");
#endif
	Long Nr = a.nrows(), Nc = a.ncols(), k1 = 0;
	for (Long j = 0; j < Nc; ++j) {
		Long k = j;
		for (Long i = 0; i < Nr; ++i) {
			a[k] = a1[k1];
			k += Nr; ++k1;
		}
	}
}

template <class T, class T1, SLS_IF(
	is_Dmat<T>() && is_Dmat<T1>())>
inline void copy(T &a, const T1 &a1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(a, a1))
		SLS_ERR("wrong shape!");
#endif
	T *p = a.ptr(), *p1 = a1.ptr();
	Long lda = a.lda(), lda1 = a1.lda();
	Long Nr = a.nrows(), Nc = a.ncols();
	for (Long j = 0; j < Nr; ++j) {
		veccpy(p, p1, Nc);
		p += lda; p1 += lda1;
	}
}

template <class T, class T1, SLS_IF(
	is_Dcmat<T>() && is_Dcmat<T1>())>
inline void copy(T &a, const T1 &a1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(a, a1))
		SLS_ERR("wrong shape!");
#endif
	T *p = a.ptr(), *p1 = a1.ptr();
	Long lda = a.lda(), lda1 = a1.lda();
	Long Nr = a.nrows(), Nc = a.ncols();
	for (Long j = 0; j < Nc; ++j) {
		veccpy(p, p1, Nr);
		p += lda; p1 += lda1;
	}
}

template <class T, class T1, SLS_IF(
	is_Dvector<T>() && is_Dvector<T1>() &&
	is_promo<contain_type<T>, contain_type<T1>>())>
inline void copy(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	T *p_end = &v.end();
	T1 *p1 = v1.ptr();
	Long step = v.step(), step1 = v1.step();
	for (T *p = v.ptr(); p <= p_end; p += step) {
		*p = *p1;
		p1 += step1;
	}
}

template <class T, class T1, SLS_IF(
	is_dense_vec<T>() && is_Dvector<T1>() &&
	is_promo<contain_type<T>, contain_type<T1>>())>
inline void copy(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	auto *p_end = &v.end();
	auto *p1 = v1.ptr();
	Long step1 = v1.step();
	for (auto p = v.ptr(); p <= p_end; ++p) {
		*p = *p1;
		p1 += step1;
	}
}

template <class T, class T1, SLS_IF(
	is_Dvector<T>() && is_dense_vec<T1>() &&
	is_promo<contain_type<T>, contain_type<T1>>())>
inline void copy(T &v, const T1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (!shape_cmp(v, v1))
		SLS_ERR("wrong shape!");
#endif
	T *p1_end = &v1.end();
	T1 *p = v.ptr();
	Long step = v.step();
	for (T *p1 = v1.ptr(); p1 <= p1_end; ++p1) {
		*p = *p1;
		p += step;
	}
}

} // namespace slisc
