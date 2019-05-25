// copy data from one container to another
// includes container shape checking
// all `container = container` should be implemented using copy

#include "meta.h"

namespace slisc {
// memory set and copy
template<class T>
inline void vecset(T *dest, const T &val, Long_I n)
{
	for (Long i = 0; i < n; ++i)
		dest[i] = val;
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
	major<T>() == major<T1>())>
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
	major<T>() == 'c' && major<T1>() == 'r')>
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
	major<T>() == 'r' && major<T1>() == 'c')>
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

} // namespace slisc
