// row-major 3D container
#pragma once
#include "vector.h"

namespace slisc {
// 3D Matrix Class

template <class T>
class Mat3d : public Vbase<T>
{
private:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_N1;
    Long m_N2;
    Long m_N3;
    Mat3d(); // default constructor: uninitialized
public:
    using Base::operator();
    using Base::ptr;
    using Base::operator=;
    Mat3d(Long_I N1, Long_I N2, Long_I N3);
    Mat3d(Long_I N1, Long_I N2, Long_I N3, const T &a);
    Mat3d(const Mat3d &rhs);   // Copy constructor
    Mat3d & operator=(const Mat3d &rhs);    // copy assignment
    template <class Tmat3, SLS_IF(is_dense_mat3<Tmat3>())>
    Mat3d & operator=(const Tmat3 &rhs);
#ifdef _CUSLISC_
    Mat3d & operator=(const Gmat3d<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
    void operator<<(Mat3d &rhs); // move data and rhs.resize(0, 0, 0)
    void resize(Long_I N1, Long_I N2, Long_I N3);
    template <class T1>
    void resize(const Mat3d<T1> &a);
    T & operator()(Long_I i, Long_I j, Long_I k);    //subscripting: pointer to row i
    const T & operator()(Long_I i, Long_I j, Long_I k) const;
    Long n1() const;
    Long n2() const;
    Long n3() const;
};

template <class T>
Mat3d<T>::Mat3d() {}

template <class T>
Mat3d<T>::Mat3d(Long_I N1, Long_I N2, Long_I N3) : Base(N1*N2*N3), m_N1(N1), m_N2(N2), m_N3(N3) {}

template <class T>
Mat3d<T>::Mat3d(Long_I N1, Long_I N2, Long_I N3, const T &s) : Mat3d(N1, N2, N3)
{ *this = s; }

template <class T>
Mat3d<T>::Mat3d(const Mat3d<T> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Mat3d<T> & Mat3d<T>::operator=(const Mat3d<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T> template <class Tmat3, SLS_IF0(is_dense_mat3<Tmat3>())>
inline Mat3d<T> & Mat3d<T>::operator=(const Tmat3 &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline void Mat3d<T>::operator<<(Mat3d<T> &rhs)
{
    m_N1 = rhs.m_N1; m_N2 = rhs.m_N2; m_N3 = rhs.m_N3;
    rhs.m_N1 = rhs.m_N2 = rhs.m_N3 = 0;
    Base::operator<<(rhs);
}

template <class T>
inline void Mat3d<T>::resize(Long_I N1, Long_I N2, Long_I N3)
{
    if (N1 != m_N1 || N2 != m_N2 || N3 != m_N3) {
        Base::resize(N1*N2*N3);
        m_N1 = N1; m_N2 = N2; m_N3 = N3;
    }
}

template <class T>
template <class T1>
inline void Mat3d<T>::resize(const Mat3d<T1> &a) { resize(a.n1(), a.n2(), a.n3()); }

template <class T>
inline T & Mat3d<T>::operator()(Long_I i, Long_I j, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 || k < 0 || k >= m_N3)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_N2*m_N3*i + m_N3*j + k];
}

template <class T>
inline const T & Mat3d<T>::operator()(Long_I i, Long_I j, Long_I k) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 || k < 0 || k >= m_N3)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_N2*m_N3*i + m_N3*j + k];
}

template <class T>
inline Long Mat3d<T>::n1() const { return m_N1; }

template <class T>
inline Long Mat3d<T>::n2() const { return m_N2; }

template <class T>
inline Long Mat3d<T>::n3() const { return m_N3; }
} // namespace slisc
