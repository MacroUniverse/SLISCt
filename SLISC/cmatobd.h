// overlapping block diagonal matrix (overlap by one element)
// such as the kinetic matrix for 1D FEDVR grid
// assume it is hamiltonian
#pragma once
#include "scmat.h"

namespace slisc {

template <class T>
class CmatObd
{
protected:
	Cmat3d<T> m_data;
	Long m_N1; // m_N2 = m_N2 = (blk_size - 1) * Nblk + 1
	CmatObd() {}; // default constructor forbidden
public:
	CmatObd(Long_I blk_size, Long_I Nblk);
	CmatObd(const MatCoo<T> &a, Long_I blk_size, Long_I Nblk);
	CmatObd(const Cmat3d<T> &a);
	const T * ptr() const;
	T * ptr();
	Long n1() const;
	Long n2() const;
	Long size() const;
	Long n0() const; // n0() = m_data.n1() = m_data.n2()
	Long nblk() const; // m_data.n3()
	const T operator()(Long_I i, Long_I j) const;
};

template <class T>
CmatObd<T>::CmatObd(Long_I blk_size, Long_I Nblk)
	: m_data(blk_size, blk_size, Nblk), m_N1((blk_size - 1) * Nblk + 1)
{}

// convert Mcoo matrix to MatOdb matrix
template <class T>
CmatObd<T>::CmatObd(const MatCoo<T> &a, Long_I blk_size, Long_I Nblk)
	: CmatObd(blk_size, Nblk)
{
#ifdef SLS_CHECK_SHAPE
	if (a.n1() != m_N1 || a.n2() != m_N1)
		SLS_ERR("wrong shape!");
#endif
	for (Long k = 0; k < a.nnz(); ++k) {
		Long i = a.row(k), j = a.col(k);
		Long N = n0() - 1, Nblk = m_data.n3();
		Long iblk = i / N, jblk = j / N;
		Long m = i % N;
		if (iblk == jblk) {
			if (iblk == Nblk)
				m_data(N, N, Nblk - 1) = a[k];
			else if (i == j && m == 0 && iblk > 0)
				 m_data(0, 0, iblk) =
					 m_data(N, N, iblk - 1) = T(a[k] * 0.5);
			else
				m_data(m, j % N, iblk) = a[k];
			continue;
		}
		else if (jblk == iblk - 1) {
			if (m == 0) {
				m_data(N, j % N, jblk) = a[k];
				continue;
			}
		}
		else if (jblk == iblk + 1) {
			Long n = j % N;
			if (n == 0) {
				m_data(m, N, iblk) = a[k];
				continue;
			}
		}
		SLS_ERR("element in a out of block!");
	}
}

template<class T>
inline CmatObd<T>::CmatObd(const Cmat3d<T>& a)
	: CmatObd(a.n1(), a.n3())
{
#ifdef SLS_CHECK_SHAPE
	if (a.n1() != a.n2())
		SLS_ERR("input must be a square matrix!");
	Long end = n0() - 1;
	for (Long k = a.n3() - 2; k >= 0; --k) {
		if (a(end, end, k) != a(0, 0, k+1))
			SLS_ERR("overlapping elements not equal!");
	}
#endif
	m_data = a;
}

template<class T>
const T * CmatObd<T>::ptr() const
{
	return m_data.ptr();
}

template<class T>
T * CmatObd<T>::ptr()
{
	return m_data.ptr();
}

template <class T>
Long CmatObd<T>::n1() const
{
	return m_N1;
}

template <class T>
Long CmatObd<T>::n2() const
{
	return m_N1;
}

template <class T>
Long CmatObd<T>::size() const
{
	return m_N1 * m_N1;
}

template<class T>
inline Long CmatObd<T>::n0() const
{
	return m_data.n1();
}

template<class T>
inline Long CmatObd<T>::nblk() const
{
	return m_data.n3();
}

template<class T>
const T CmatObd<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
	if (i < 0 || i >= m_N1 || j < 0 || j >= m_N1)
		SLS_ERR("out of bound!");
#endif
	Long N = n0() - 1, Nblk = m_data.n3();
	Long iblk = i / N, jblk = j / N;
	Long m = i % N;
	if (iblk == jblk) {
		if (iblk == Nblk)
			return m_data(N, N, Nblk - 1);
		else if (i == j && m == 0 && iblk > 0)
			return 2 * m_data(0, 0, iblk);
		return m_data(m, j % N, iblk);
	}
	else if (jblk == iblk - 1) {
		if (m == 0)
			return m_data(N, j % N, jblk);
	}
	else if (jblk == iblk + 1) {
		Long n = j % N;
		if (n == 0)
			return m_data(m, N, iblk);
	}
	return 0;
}

} // namespace slisc
