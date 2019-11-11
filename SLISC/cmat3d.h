// column-major 3D container
#pragma once
#include "vector.h"
#include "scalar_arith.h"

namespace slisc {
// 3D Matrix Class

template <class T>
class Cmat3d : public Vbase<T>
{
private:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_N1;
    Long m_N2;
    Long m_N3;
    Cmat3d(); // default constructor: uninitialized
public:
    typedef T value_type;
    using Base::operator();
    using Base::ptr;
    using Base::operator=;
    Cmat3d(Long_I N1, Long_I N2, Long_I N3);
    Cmat3d(Long_I N1, Long_I N2, Long_I N3, const T &a);
    Cmat3d(const Cmat3d &rhs);   // Copy constructor
    Cmat3d & operator=(const Cmat3d &rhs);    // copy assignment
    template <class T1>
    Cmat3d & operator=(const Cmat3d<T1> &rhs);
#ifdef _CUSLISC_
    Cmat3d & operator=(const Gmat3d<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
#ifdef _CUSLISC_
    Cmat3d & operator=(const Gcmat3d<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
    void operator<<(Cmat3d &rhs); // move data and rhs.resize(0, 0, 0)
    void resize(Long_I N1, Long_I N2, Long_I N3);
    template <class T1>
    void resize(const Cmat3d<T1> &a);
    void resize_cpy(Long_I N1, Long_I N2, Long_I N3);
    T & operator()(Long_I i, Long_I j, Long_I k);    //subscripting: pointer to row i
    const T & operator()(Long_I i, Long_I j, Long_I k) const;
    Long n1() const;
    Long n2() const;
    Long n3() const;
};

template <class T>
Cmat3d<T>::Cmat3d() {}

template <class T>
Cmat3d<T>::Cmat3d(Long_I N1, Long_I N2, Long_I N3) : Base(N1*N2*N3), m_N1(N1), m_N2(N2), m_N3(N3) {}

template <class T>
Cmat3d<T>::Cmat3d(Long_I N1, Long_I N2, Long_I N3, const T &s) : Cmat3d(N1, N2, N3)
{ *this = s; }

template <class T>
Cmat3d<T>::Cmat3d(const Cmat3d<T> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Cmat3d<T> & Cmat3d<T>::operator=(const Cmat3d<T> &rhs)
{
    return operator=<T>(rhs);
}

template <class T> template <class T1>
inline Cmat3d<T> & Cmat3d<T>::operator=(const Cmat3d<T1> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline void Cmat3d<T>::operator<<(Cmat3d<T> &rhs)
{
    m_N1 = rhs.m_N1; m_N2 = rhs.m_N2; m_N3 = rhs.m_N3;
    rhs.m_N1 = rhs.m_N2 = rhs.m_N3 = 0;
    Base::operator<<(rhs);
}

template <class T>
inline void Cmat3d<T>::resize(Long_I N1, Long_I N2, Long_I N3)
{
    if (N1 != m_N1 || N2 != m_N2 || N3 != m_N3) {
        Base::resize(N1*N2*N3);
        m_N1 = N1; m_N2 = N2; m_N3 = N3;
    }
}

template <class T>
template <class T1>
inline void Cmat3d<T>::resize(const Cmat3d<T1> &a) { resize(a.n1(), a.n2(), a.n3()); }

template<class T>
inline void Cmat3d<T>::resize_cpy(Long_I N1, Long_I N2, Long_I N3)
{
    // todo: not optimized
    if (N1 != m_N1 || N2 != m_N2 || N3 != m_N3) {
        Long min_N1 = min(m_N1, N1), min_N2 = min(m_N2, N2), min_N3 = min(m_N3, N3);
        Cmat3d<T> new_mat(N1, N2, N3, 0);
        slice(new_mat, 0, min_N1, 0, min_N2, 0, min_N3) =
            slice(*this, 0, min_N1, 0, min_N2, 0, min_N3);
        *this << new_mat;
    }
}

template <class T>
inline T & Cmat3d<T>::operator()(Long_I i, Long_I j, Long_I k)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 || k < 0 || k >= m_N3)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i + m_N1*j + m_N1*m_N2*k];
}

template <class T>
inline const T & Cmat3d<T>::operator()(Long_I i, Long_I j, Long_I k) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 || k < 0 || k >= m_N3)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i + m_N1*j + m_N1*m_N2*k];
}

template <class T>
inline Long Cmat3d<T>::n1() const {
    return m_N1;
}

template <class T>
inline Long Cmat3d<T>::n2() const {
    return m_N2;
}

template <class T>
inline Long Cmat3d<T>::n3() const {
    return m_N3;
}
} // namespace slisc
