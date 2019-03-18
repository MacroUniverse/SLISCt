#pragma once
#include "vector.h"

namespace slisc {
// === slice to Vector ===

// #### IMPORTANT ####
// the sliced vector "v" must not be deleted or resized
// "v" must be empty or a sliced vector before using any function below
// use "slice_reset(v)" before "v" is destroyed

// contiguous slice, set size
template <class T>
inline void slice(Vector<T> &v, const T *ptr, Long_I N)
{
	v.m_p = (T *)ptr; v.m_N = N;
}
	
// contiguous slice, size unchanged
template <class T>
inline void slice(Vector<T> &v, const T *ptr)
{
	v.m_p = (T *)ptr;
}

// set size (m_N) of the vector, pointer unchanged
template <class T>
inline void slice_resize(Vector<T> &v, Long_I N)
{
#ifdef SLS_CHECK_SHAPE
	if (N <= 0) error("illegal N!");
#endif
	v.m_N = N;
}

// reset vector to normal, with zero size
template <class T>
inline void slice_reset(Vector<T> &v)
{
	v.m_p = nullptr; v.m_N = 0;
}
}
