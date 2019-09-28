// overlapping block diagonal matrix (overlap by one element)
// such as the kinetic matrix for 1D FEDVR grid
// first block and last block has one less element in each dimension

#pragma once
#include "scmat.h"

namespace slisc {

template <class T>
class CmatObd
{
protected:
    Cmat3d<T> m_data;
    Long m_N1; // m_N2 = m_N2 = (blk_size - 1) * Nblk - 1
    CmatObd() {}; // default constructor forbidden
public:
    typedef T value_type;
    CmatObd(Long_I blk_size, Long_I Nblk);
    const T &operator()(Long_I i) const; // m_data[i]
    T &operator()(Long_I i);
    Long find(Long_I i, Long_I j);
    template <class T1, SLS_IF(is_promo<T, T1>())>
    CmatObd &operator=(const CmatObd<T1> &rhs);
    template <class T1, SLS_IF(is_promo<T, T1>())>
    CmatObd &operator=(const MatCoo<T1> &rhs);
    template <class T1, SLS_IF(is_promo<T, T1>())>
    CmatObd &operator=(const Cmat3d<T1> &a);
    const T * ptr() const; // not the first element!
    T * ptr();
    Long n1() const;
    Long n2() const;
    Long size() const;
    Long nnz() const;
    const Cmat3d<T> &cmat3() const;
    Cmat3d<T> &cmat3();
    Long n0() const; // n0() = m_data.n1() = m_data.n2()
    Long nblk() const; // m_data.n3()
    const T operator()(Long_I i, Long_I j) const;
    void resize(Long_I blk_size, Long_I Nblk);
};

template <class T>
CmatObd<T>::CmatObd(Long_I blk_size, Long_I Nblk)
    : m_data(blk_size, blk_size, Nblk), m_N1((blk_size - 1) * Nblk - 1)
{
    Long step = SQR(n0());
    // set the first overlapped element to 0
    if (m_data.size() > 0)
        vecset(m_data.ptr() + step - 1, 0, Nblk - 1, step);
}

template<class T>
const T & CmatObd<T>::operator()(Long_I i) const
{
    return m_data[i];
}

template<class T>
T & CmatObd<T>::operator()(Long_I i)
{
    return m_data[i];
}

template<class T>
inline Long CmatObd<T>::find(Long_I i1, Long_I i2)
{
#ifdef SLS_CHECK_BOUNDS
    if (i1 < 0 || i1 >= m_N1 || i2 < 0 || i2 >= m_N1)
        SLS_ERR("out of bound!");
#endif
    Long i = i1 + 1; Long j = i2 + 1;
    Long N = n0() - 1, Nblk = m_data.n3();
    Long iblk = i / N, jblk = j / N;
    Long m = i % N;
    Long step2 = n0(), step3 = step2 * step2;
    if (iblk == jblk) {
        if (iblk == Nblk)
            return N + step2 * N + step3 * (Nblk - 1);
        else if (i == j && m == 0 && iblk > 0)
            return step3 * iblk;
        return m + step2 * (j % N) + step3 * iblk;
    }
    else if (jblk == iblk - 1) {
        if (m == 0)
            return N + step2 * (j % N) + step3 * jblk;
    }
    else if (jblk == iblk + 1) {
        Long n = j % N;
        if (n == 0)
            return m + step2 * N + step3 * iblk;
    }
    SLS_ERR("element out of block!");
    return -1;
}

template <class T>
template <class T1, SLS_IF0(is_promo<T, T1>())>
CmatObd<T> &CmatObd<T>::operator=(const CmatObd<T1> &a)
{
    m_data = a.cmat3();
    return *this;
}

// convert Mcoo matrix to MatOdb matrix
template <class T>
template <class T1, SLS_IF0(is_promo<T, T1>())>
CmatObd<T> &CmatObd<T>::operator=(const MatCoo<T1> &a)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(*this, a))
        SLS_ERR("wrong shape!");
#endif
    m_data = 0;
    for (Long k = 0; k < a.nnz(); ++k) {
        Long i = a.row(k) + 1, j = a.col(k) + 1;
        Long N = n0() - 1, Nblk = m_data.n3();
        Long iblk = i / N, jblk = j / N;
        Long m = i % N;
        if (iblk == jblk) {
            if (iblk == Nblk)
                m_data(N, N, Nblk - 1) = a[k];
            else if (i == j && m == 0 && iblk > 0)
                 m_data(0, 0, iblk) = a[k];
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
        SLS_ERR("element out of block!");
    }
    return *this;
}

template<class T>
template <class T1, SLS_IF0(is_promo<T, T1>())>
CmatObd<T> &CmatObd<T>::operator=(const Cmat3d<T1>& rhs)
{
    m_data = rhs;
    // set the first overlapped element to 0
    Long step = SQR(n0());
    vecset(m_data.ptr() + step - 1, 0, nblk() - 1, step);
    return *this;
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
Long CmatObd<T>::nnz() const
{
    Long N0 = n0(), Nblk = nblk();
    return (N0*N0 - 1)*Nblk - 4 * N0 + 3;
}

template<class T>
const Cmat3d<T>& CmatObd<T>::cmat3() const
{
    return m_data;
}

template<class T>
Cmat3d<T>& CmatObd<T>::cmat3()
{
    return m_data;
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
const T CmatObd<T>::operator()(Long_I i1, Long_I i2) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i1 < 0 || i1 >= m_N1 || i2 < 0 || i2 >= m_N1)
        SLS_ERR("out of bound!");
#endif
    Long i = i1 + 1, j = i2 + 1;
    Long N = n0() - 1, Nblk = m_data.n3();
    Long iblk = i / N, jblk = j / N;
    Long m = i % N;
    if (iblk == jblk) {
        if (iblk == Nblk)
            return m_data(N, N, Nblk - 1);
        else if (i == j && m == 0 && iblk > 0)
            return m_data(0, 0, iblk);
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
    return T(0);
}

template<class T>
inline void CmatObd<T>::resize(Long_I blk_size, Long_I Nblk)
{
    m_data.resize(blk_size, blk_size, Nblk);
    m_N1 = (blk_size - 1) * Nblk - 1;
    Long step = SQR(n0());
    // set the first overlapped element to 0
    vecset(m_data.ptr() + step - 1, 0, Nblk - 1, step);
}

} // namespace slisc
