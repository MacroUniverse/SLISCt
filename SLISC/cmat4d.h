// column-major 4D container
#pragma once
#include "vector.h"

namespace slisc {
// 3D Matrix Class

template <class T>
class Cmat4d : public Vbase<T>
{
private:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_N1, m_N2, m_N3, m_N4;
    Cmat4d(); // default constructor: uninitialized
public:
    typedef T value_type;
    using Base::operator();
    using Base::ptr;
    using Base::operator=;
    Cmat4d(Long_I N1, Long_I N2, Long_I N3, Long_I N4);
    Cmat4d(Long_I N1, Long_I N2, Long_I N3, Long_I N4, const T &a);
    Cmat4d(const Cmat4d &rhs);   // Copy constructor
    Cmat4d & operator=(const Cmat4d &rhs);    // copy assignment
    template <class T1>
    Cmat4d & operator=(const Cmat4d<T1> &rhs);
#ifdef _CUSLISC_
    Cmat4d & operator=(const Gmat3d<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
    void operator<<(Cmat4d &rhs); // move data and rhs.resize(0, 0, 0)
    void resize(Long_I N1, Long_I N2, Long_I N3, Long_I N4);
    template <class T1>
    void resize(const Cmat4d<T1> &a);
    T & operator()(Long_I i, Long_I j, Long_I k, Long_I l);
    const T & operator()(Long_I i, Long_I j, Long_I k, Long_I l) const;
    Long n1() const;
    Long n2() const;
    Long n3() const;
    Long n4() const;
};

template <class T>
Cmat4d<T>::Cmat4d() {}

template <class T>
Cmat4d<T>::Cmat4d(Long_I N1, Long_I N2, Long_I N3, Long_I N4)
    : Base(N1*N2*N3*N4), m_N1(N1), m_N2(N2), m_N3(N3), m_N4(N4) {}

template <class T>
Cmat4d<T>::Cmat4d(Long_I N1, Long_I N2, Long_I N3, Long_I N4, const T &s)
    : Cmat4d(N1, N2, N3, N4) {
    *this = s;
}

template <class T>
Cmat4d<T>::Cmat4d(const Cmat4d<T> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Cmat4d<T> & Cmat4d<T>::operator=(const Cmat4d<T> &rhs)
{
    return operator=<T>(rhs);
}

template <class T> template <class T1>
inline Cmat4d<T> & Cmat4d<T>::operator=(const Cmat4d<T1> &rhs)
{
#ifdef SLS_CHECK_SHAPE
    if (m_N1 != rhs.n1() || m_N2 != rhs.n2() || m_N3 != rhs.n3())
        SLS_ERR("wrong shape!");
#endif
    copy(*this, rhs);
    return *this;
}

template <class T>
inline void Cmat4d<T>::operator<<(Cmat4d<T> &rhs)
{
    m_N1 = rhs.m_N1; m_N2 = rhs.m_N2; m_N3 = rhs.m_N3;
    rhs.m_N1 = rhs.m_N2 = rhs.m_N3 = 0;
    Base::operator<<(rhs);
}

template <class T>
inline void Cmat4d<T>::resize(Long_I N1, Long_I N2, Long_I N3, Long_I N4)
{
    if (N1 != m_N1 || N2 != m_N2 || N3 != m_N3 || N4 != m_N4) {
        Base::resize(N1*N2*N3*N4);
        m_N1 = N1; m_N2 = N2; m_N3 = N3; m_N4 = N4;
    }
}

template <class T>
template <class T1>
inline void Cmat4d<T>::resize(const Cmat4d<T1> &a) { resize(a.n1(), a.n2(), a.n3()); }

template <class T>
inline T & Cmat4d<T>::operator()(Long_I i, Long_I j, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 ||
        k < 0 || k >= m_N3 || l < 0 || l >= m_N4)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    Long N1N2 = m_N1 * m_N2;
    return m_p[i + m_N1*j + N1N2 *k + N1N2 *m_N3*l];
}

template <class T>
inline const T & Cmat4d<T>::operator()(Long_I i, Long_I j, Long_I k, Long_I l) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 ||
        k < 0 || k >= m_N3 || l < 0 || l >= m_N4)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    Long N1N2 = m_N1 * m_N2;
    return m_p[i + m_N1 * j + N1N2 * k + N1N2 * m_N3*l];
}

template <class T>
inline Long Cmat4d<T>::n1() const {
    return m_N1;
}

template <class T>
inline Long Cmat4d<T>::n2() const {
    return m_N2;
}

template <class T>
inline Long Cmat4d<T>::n3() const {
    return m_N3;
}

template <class T>
inline Long Cmat4d<T>::n4() const {
    return m_N4;
}
} // namespace slisc
