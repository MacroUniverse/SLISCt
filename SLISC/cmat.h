// column-major matrix container
#pragma once
#include "matrix.h"

namespace slisc {

// Column major Matrix Class
// don't derive from Matrix<T>

template <class T>
class Cmat : public Vbase<T>
{
protected:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_Nr, m_Nc;
    Cmat(); // default constructor: uninitialized
public:
    using Base::ptr;
    using Base::operator();
    using Base::operator=;
    Cmat(Long_I Nr, Long_I Nc);
    Cmat(Long_I Nr, Long_I Nc, const T &s);    //Initialize to constant
    Cmat(Long_I Nr, Long_I Nc, const T *ptr);    // Initialize to array
    Cmat(const Cmat &rhs);        // Copy constructor
    Cmat & operator=(const Cmat &rhs);    // copy assignment
    template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
    Cmat & operator=(const Tmat &rhs);
    Cmat & operator=(const T &rhs);
    template <class T1>
    Cmat & operator=(const MatCoo<T1> &rhs);
    template <class T1>
    Cmat & operator=(const MatCooH<T1> &rhs);
#ifdef _CUSLISC_
    Cmat & operator=(const Gcmat<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
    void operator<<(Cmat &rhs); // move data and rhs.resize(0, 0)
    T& operator()(Long_I i, Long_I j);    // double indexing
    const T& operator()(Long_I i, Long_I j) const;
    Long n1() const;
    Long n2() const;
    void resize(Long_I Nr, Long_I Nc); // resize (contents not preserved)
    template <class T1>
    void resize(const Cmat<T1> &a);
    void resize_cpy(Long_I N1, Long_I N2);
};

template <class T>
Cmat<T>::Cmat() {}

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc) : Base(Nr*Nc), m_Nr(Nr), m_Nc(Nc) {}

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc, const T &s) : Cmat(Nr, Nc)
{ *this = s; }

template <class T>
Cmat<T>::Cmat(Long_I Nr, Long_I Nc, const T *ptr) : Cmat(Nr, Nc)
{ memcpy(m_p, ptr, m_N*sizeof(T)); }

template <class T>
Cmat<T>::Cmat(const Cmat<T> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, "
        "use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Cmat<T> & Cmat<T>::operator=(const Cmat<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T> template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
inline Cmat<T> & Cmat<T>::operator=(const Tmat &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline Cmat<T> & Cmat<T>::operator=(const T &rhs)
{
    vecset(m_p, rhs, m_N);
    return *this;
}

template <class T> template <class T1>
inline Cmat<T> & Cmat<T>::operator=(const MatCoo<T1> &rhs)
{
    return coo2dense(*this, rhs);
}

template <class T> template <class T1>
inline Cmat<T> & Cmat<T>::operator=(const MatCooH<T1> &rhs)
{
    return cooh2dense(*this, rhs);
}

template <class T>
inline void Cmat<T>::operator<<(Cmat<T> &rhs)
{
    m_Nr = rhs.m_Nr; m_Nc = rhs.m_Nc;
    rhs.m_Nr = rhs.m_Nc = 0;
    Base::operator<<(rhs);
}

template <class T>
inline T & Cmat<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i+m_Nr*j];
}

template <class T>
inline const T & Cmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i+m_Nr*j];
}

template <class T>
inline Long Cmat<T>::n1() const
{ return m_Nr; }

template <class T>
inline Long Cmat<T>::n2() const
{ return m_Nc; }

template <class T>
inline void Cmat<T>::resize(Long_I Nr, Long_I Nc)
{
    if (Nr != m_Nr || Nc != m_Nc) {
        Base::resize(Nr*Nc);
        m_Nr = Nr; m_Nc = Nc;
    }
}

template<class T>
inline void Cmat<T>::resize_cpy(Long_I N1, Long_I N2)
{
    Long N = N1 * N2;
    if (N1 != m_Nr || N2 != m_Nc) {
        if (m_N == 0 || N == 0) {
            resize(N1, N2);
        }
        else {
            T *old_p = m_p;
            m_p = new T[N];
            if (N > m_N) {
                for (Long j = 0; j < m_Nc; ++j) {
                    veccpy(m_p + N1*j, old_p + m_Nr*j, m_Nr);
                    vecset(m_p + N1*j + m_Nr, 0, N1 - m_Nr);
                }
                vecset(m_p + N1 * m_Nc, 0, N1 * (N2 - m_Nc));
            }
            else {// N < m_N
                for (Long j = 0; j < N2; ++j) {
                    veccpy(m_p + N1 * j, old_p + m_Nr * j, N1);
                }
            }
            m_Nr = N1; m_Nc = N2; m_N = N;
            delete[] old_p;
        }
    }
}

template <class T>
template <class T1>
inline void Cmat<T>::resize(const Cmat<T1> &a)
{ resize(a.n1(), a.n2()); }

}
