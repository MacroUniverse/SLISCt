#include "vector.h"

namespace slisc {
	// Vector Class
	// the sliced vector cannot be deleted or resized, use at your own risk
	// template <class T, SLS_IF(is_scalar<T>)>
	template <class T>
	inline void slice(Vector<T> &v, const T *ptr)
	{
		v.m_p = (T *)ptr;
	}

	template <class T>
	inline void slice(Vector<T> &v, const T *ptr, Long_I N)
	{
		v.m_p = (T *)ptr; v.m_N = N;
	}

	template <class T>
	inline void rm_slice(Vector<T> &v)
	{
		v.m_p = nullptr; v.m_N = 0;
	}
}
