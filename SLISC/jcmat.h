// matrix slicing with both `step1` and `step2`
// Jcmat can be used for both col and row major
#pragma once
#include "vector.h"

namespace slisc {

template <class T>
class Jcmat
{
private:
    T *m_p;
    Long m_N;
    Long m_N1, m_N2;
    Long m_step1, m_step2; // a(i,j) = m_p + m_step1*i + m_step2*j
public:
    typedef T value_type;
    Jcmat();
    Jcmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I rstep, Long_I cstep);
    void set(const T *ptr, Long_I Nr, Long_I Nc, Long_I rstep, Long_I cstep);

    // === Cmat member functions ===
    Jcmat & operator=(const Jcmat &rhs);    // copy assignment
    template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
    Jcmat & operator=(const Tmat &rhs);
    Jcmat & operator=(const T &rhs);
    //template <class T1>
    //Cmat & operator=(const MatCoo<T1> &rhs);
    //template <class T1>
    //Cmat & operator=(const MatCooH<T1> &rhs);
    T& operator[](Long_I i);    // single indexing (inefficient)
    const T& operator[](Long_I i) const;
    T& operator()(Long_I i);    // same as operator[]
    const T& operator()(Long_I i) const;
    T& operator()(Long_I i, Long_I j);    // double indexing
    const T& operator()(Long_I i, Long_I j) const;
    Long n1() const;
    Long n2() const;
    Long step1() const;
    Long step2() const;
    Long size() const;

    const T *ptr() const;
    T *ptr();
};

template <class T>
Jcmat<T>::Jcmat() {}

template <class T>
Jcmat<T>::Jcmat(const T *ptr, Long_I Nr, Long_I Nc, Long_I step1, Long_I step2)
    : m_p((T *)ptr), m_N1(Nr), m_N2(Nc), m_N(Nr*Nc), m_step1(step1), m_step2(step2)
{}

template <class T>
void Jcmat<T>::set(const T *ptr, Long_I Nr, Long_I Nc, Long_I step1, Long_I step2)
{
    m_p = (T *)ptr; m_N1 = Nr; m_N2 = Nc;
    m_step1 = step1; m_step2 = step2;
    m_N = Nr * Nc;
}

template <class T>
inline Jcmat<T> & Jcmat<T>::operator=(const Jcmat<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
Jcmat<T> & Jcmat<T>::operator=(const Tmat &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline Jcmat<T> & Jcmat<T>::operator=(const T &rhs)
{
    T *p = m_p;
    for (Long j = 0; j < m_N2; ++j) {
        for (Long i = 0; i < m_N1; ++i) {
            *p = rhs;
            p += m_step1;
        }
        p += m_step2;
    }
    return *this;
}

template <class T>
inline T & Jcmat<T>::operator[](Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * (i % m_N1) + m_step2 * (i / m_N1)];
}

template <class T>
inline const T & Jcmat<T>::operator[](Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * (i % m_N1) + m_step2 * (i / m_N1)];
}

template <class T>
inline T & Jcmat<T>::operator()(Long_I i)
{
    return operator[](i);
}

template <class T>
inline const T & Jcmat<T>::operator()(Long_I i) const
{
    return operator[](i);
}

template <class T>
inline T & Jcmat<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * i + m_step2 * j];
}

template <class T>
inline const T & Jcmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * i + m_step2 * j];
}

template <class T>
inline Long Jcmat<T>::n1() const
{
    return m_N1;
}

template <class T>
inline Long Jcmat<T>::n2() const
{
    return m_N2;
}

template <class T>
inline Long Jcmat<T>::step1() const
{
    return m_step1;
}

template <class T>
inline Long Jcmat<T>::step2() const
{
    return m_step2;
}

template<class T>
inline Long Jcmat<T>::size() const
{
    return m_N;
}

template<class T>
inline const T * Jcmat<T>::ptr() const
{
    return m_p;
}

template<class T>
inline T * Jcmat<T>::ptr()
{
    return m_p;
}

} // namespace slisc
