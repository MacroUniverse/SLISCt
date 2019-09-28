#pragma once
#include "vector.h"

namespace slisc {

// contiguous slice matrix class (column major)
// const Scmat is only top level const
// use Scmat_c for low level const
template <class T>
class Scmat : public Svector<T>
{
public:
    typedef Svector<T> Base;
    typedef T value_type;
    using Base::m_p;
    using Base::m_N;
    Long m_Nr, m_Nc;
    Scmat();
    Scmat(Long_I Nr, Long_I Nc);
    Scmat(T *ptr, Long_I Nr, Long_I Nc);

    // === Cmat functions ===
    Scmat & operator=(const Scmat &rhs);    // copy assignment
    template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
    Scmat & operator=(const Tmat &rhs);
    Scmat & operator=(const T &rhs);
    template <class T1>
    Scmat & operator=(const MatCoo<T1> &rhs);
    template <class T1>
    Scmat & operator=(const MatCooH<T1> &rhs);
    T& operator()(Long_I i, Long_I j) const; // double indexing
    Long n1() const;
    Long n2() const;

    // resize() is a bad idea, don't try to create it!

    // There is no upper bound checking of N, use with care
    void set_size(Long_I Nr, Long_I Nc);
    void set_ptr(T *ptr);
    void set(T *ptr, Long_I Nr, Long_I Nc);
    void next(); // m_ptr += m_N
    void last(); // m_ptr -= m_N
    void shift(Long_I N); // m_ptr += N;
    ~Scmat();
};

template <class T>
inline Scmat<T>::Scmat() {}

template <class T>
inline Scmat<T>::Scmat(Long_I Nr, Long_I Nc)
    : Base(Nr*Nc), m_Nr(Nr), m_Nc(Nc) {}

template <class T>
inline Scmat<T>::Scmat(T *ptr, Long_I Nr, Long_I Nc)
    : Scmat(Nr, Nc)
{
    m_p = (T *)ptr;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const Scmat<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T> template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
inline Scmat<T> & Scmat<T>::operator=(const Tmat &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline Scmat<T> & Scmat<T>::operator=(const T &rhs)
{
    vecset(m_p, rhs, m_N);
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCoo<T1> &rhs)
{
    return coo2dense(*this, rhs);
}

template <class T> template <class T1>
inline Scmat<T> & Scmat<T>::operator=(const MatCooH<T1> &rhs)
{
    return cooh2dense(*this, rhs);
}

template <class T>
inline T & Scmat<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i + m_Nr * j];
}

template <class T>
inline Long Scmat<T>::n1() const
{
    return m_Nr;
}

template <class T>
inline Long Scmat<T>::n2() const
{
    return m_Nc;
}

template <class T>
inline void Scmat<T>::set_size(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
    if (Nr <= 0 || Nc <= 0) SLS_ERR("illegal Nr or Nc!");
#endif
    m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat<T>::set_ptr(T * ptr)
{
    m_p = (T *)ptr;
}

template <class T>
inline void Scmat<T>::set(T * ptr, Long_I Nr, Long_I Nc)
{
    m_p = (T *)ptr;
    m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat<T>::next()
{
    m_p += m_N;
}

template <class T>
inline void Scmat<T>::last()
{
    m_p -= m_N;
}

template <class T>
inline void Scmat<T>::shift(Long_I N)
{
    m_p += N;
}

template <class T>
inline Scmat<T>::~Scmat() {}


// ========== low level const of Scmat =========
template <class T>
class Scmat_c : public Svector<T>
{
public:
    typedef Svector<T> Base;
    typedef T value_type;
    using Base::m_p;
    using Base::m_N;
    Long m_Nr, m_Nc;
    Scmat_c();
    Scmat_c(Long_I Nr, Long_I Nc);
    Scmat_c(const T *ptr, Long_I Nr, Long_I Nc);

    // === Cmat functions ===
    const T& operator()(Long_I i, Long_I j) const;
    Long n1() const;
    Long n2() const;

    // resize() is a bad idea, don't try to create it!

    // There is no upper bound checking of N, use with care
    void set_size(Long_I Nr, Long_I Nc);
    void set_ptr(const T *ptr);
    void set(const T *ptr, Long_I Nr, Long_I Nc);
    void next(); // m_ptr += m_N
    void last(); // m_ptr -= m_N
    void shift(Long_I N); // m_ptr += N;
    ~Scmat_c();
};

template <class T>
inline Scmat_c<T>::Scmat_c() {}

template <class T>
inline Scmat_c<T>::Scmat_c(Long_I Nr, Long_I Nc)
    : m_Nr(Nr), m_Nc(Nc), Base(Nr*Nc) {}

template <class T>
inline Scmat_c<T>::Scmat_c(const T *ptr, Long_I Nr, Long_I Nc)
    : Scmat_c(Nr, Nc)
{
    m_p = (T *)ptr;
}

template <class T>
inline const T & Scmat_c<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[i + m_Nr * j];
}

template <class T>
inline Long Scmat_c<T>::n1() const
{
    return m_Nr;
}

template <class T>
inline Long Scmat_c<T>::n2() const
{
    return m_Nc;
}

template <class T>
inline void Scmat_c<T>::set_size(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
    if (Nr <= 0 || Nc <= 0) SLS_ERR("illegal Nr or Nc!");
#endif
    m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat_c<T>::set_ptr(const T * ptr)
{
    m_p = (T *)ptr;
}

template <class T>
inline void Scmat_c<T>::set(const T * ptr, Long_I Nr, Long_I Nc)
{
    m_p = (T *)ptr;
    m_Nr = Nr; m_Nc = Nc; m_N = Nr * Nc;
}

template <class T>
inline void Scmat_c<T>::next()
{
    m_p += m_N;
}

template <class T>
inline void Scmat_c<T>::last()
{
    m_p -= m_N;
}

template <class T>
inline void Scmat_c<T>::shift(Long_I N)
{
    m_p += N;
}

template <class T>
inline Scmat_c<T>::~Scmat_c() {}
} // namespace slisc
